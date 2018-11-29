"""click_mergevcfs commands tests."""
import os
from os.path import basename, join
import subprocess
import tempfile
import shutil
import gzip
import multiprocessing

import click

from click_mergevcfs.utils import (
    add_version,
    add_PASSED_field,
    decompose_multiallelic_record,
    rename_samples_headers,
    order_vcf_columns,
)


def merge_snvs(vcf_list, out_file, working_dir):
    """For merging snvs and indels."""
    vcf_list = order_vcf_columns(vcf_list, working_dir)
    cmd = ["vcf-merge", "--collapse", "none", "--remove-duplicates"]

    for initial_vcf in vcf_list:
        vcf_filename = basename(initial_vcf)

        # decompose multiallelic records
        decomposed_vcf = join(working_dir, "decomposed_{}".format(vcf_filename))
        decompose_multiallelic_record(in_vcf=initial_vcf, out_vcf=decomposed_vcf)

        # add 'PASSED' field under INFO. ex: PASSED_caveman,PASSED_mutect
        passed_vcf = join(working_dir, "passed_{}".format(vcf_filename))
        add_PASSED_field(in_vcf=decomposed_vcf, out_vcf=passed_vcf)

        # modify column sample name to include the caller name. ex: strelka_NORMAL.
        fixed_vcf = join(working_dir, vcf_filename)
        rename_samples_headers(in_vcf=passed_vcf, out_vcf=fixed_vcf)

        # bgzip and tabix
        if not fixed_vcf.endswith(".gz"):
            fixed_vcf += ".gz"
        subprocess.check_call(["tabix", "-f", "-p", "vcf", fixed_vcf])
        cmd.append(fixed_vcf)

    temp = tempfile.NamedTemporaryFile(suffix=".tmp.vcf", dir=working_dir, delete=False)
    fout = open(temp.name, "w")
    subprocess.check_call(cmd, stdout=fout)
    fout.close()

    # vcf-merge may create multiple ALT alleles per record
    decompose_multiallelic_record(in_vcf=temp.name, out_vcf=out_file)
    # add click_mergevcfs version to the header
    add_version(out_file)


def caveman_postprocess(
    perl_path,
    flag_script,
    in_vcf,
    out_vcf,
    normal_bam,
    tumor_bam,
    bedFileLoc,
    indelBed,
    unmatchedVCFLoc,
    reference,
    flagConfig,
    flagToVcfConfig,
    annoBedLoc,
    bin_size,
    sequencing_method,
    working_dir,
):
    """Run caveman flagging on merged vcf."""

    def split_vcf(i_vcf, bin_size, temp_dir):
        """Split large vcf files into smaller files."""
        with gzip.open(i_vcf, "r") as f_in:
            lines = [l.decode() for l in f_in.readlines()]
        header = [l for l in lines if l.startswith("#")]
        records = [l for l in lines if not l.startswith("#")]
        result = []

        for x in range(0, len(records), bin_size):
            end = min(len(records), x + bin_size)
            split_vcf_file = tempfile.NamedTemporaryFile(
                prefix="split_", suffix=".vcf", dir=temp_dir, delete=False
            )
            result.append(split_vcf_file.name)
            with open(split_vcf_file.name, "w") as f_out:
                f_out.write("".join(header))
                f_out.write("".join(records[x:end]))

        # bgzip and tabix the split vcfs
        for vcf in result:
            subprocess.check_call(["bgzip", "-f", vcf])
            subprocess.check_call(["tabix", "-f", "-p", "vcf", vcf + ".gz"])
        click.echo("splitted files are: {}".format(result))
        return (result, len(range(0, len(records), bin_size)))

    def run_flagging(i_vcf, temp_dir):
        """Run caveman flagging."""
        o_vcf = tempfile.NamedTemporaryFile(
            prefix="flagged_{}".format(basename(i_vcf).strip(".vcf.gz")),
            suffix=".vcf",
            dir=temp_dir,
            delete=False,
        )

        cmd = [
            perl_path,
            flag_script,
            "-i",
            i_vcf,
            "-o",
            o_vcf.name,
            "-s",
            "HUMAN",
            "-n",
            normal_bam,
            "-m",
            tumor_bam,
            "-b",
            bedFileLoc,
            "-g",
            indelBed,
            "-umv",
            unmatchedVCFLoc,
            "-ref",
            reference + ".fai",  # Reference index (fai) from caveman help
            "-t",
            "pulldown" if sequencing_method == "TGD" else "genomic",
            "-c",
            flagConfig,
            "-v",
            flagToVcfConfig,
            "-ab",
            annoBedLoc,
            "--verbose",
        ]
        # Unicode to string
        cmd = list(map(str, cmd))
        subprocess.check_call(cmd)

        # gzip and create tbi
        subprocess.check_call(["bgzip", "-f", o_vcf.name])
        subprocess.check_call(["tabix", "-f", "-p", "vcf", o_vcf.name + ".gz"])

    def concat_vcfs(vcfs, out_vcf):
        """Concatnating a list of vcf files into one vcf."""
        result = []
        for vcf in vcfs:
            with gzip.open(vcf, "r") as f:
                lines = f.readlines()
                headers = [l.decode() for l in lines if l.startswith(b"#")]
                variants = [l.decode() for l in lines if not l.startswith(b"#")]
            result.extend(variants)

        with open(out_vcf, "w") as f:
            f.write("".join(headers))
            f.write("".join(result))

    temp_dir = tempfile.mkdtemp(
        prefix=basename(in_vcf).split(".vcf")[0], dir=working_dir
    )
    vcf_splitted_files, expected_num_file = split_vcf(in_vcf, bin_size, temp_dir)
    click.echo("split to {} vcfs, temp_dir is {}".format(expected_num_file, temp_dir))

    # Run flagging in parallel
    processes = []
    for vcf in vcf_splitted_files:
        processes.append(
            multiprocessing.Process(target=run_flagging, args=(vcf + ".gz", temp_dir))
        )
    click.echo("start multiprocessing...")
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    # Before merging, check if all split jobs exist
    split_flagged_files = [
        join(temp_dir, f)
        for f in os.listdir(temp_dir)
        if (f.startswith("flagged_split_") & f.endswith(".vcf.gz"))
    ]
    click.echo("{}, {}".format(len(split_flagged_files), expected_num_file))
    assert len(split_flagged_files) == expected_num_file

    # merge vcfs in flagged_vcfs
    click.echo("start merging...")
    unsorted_merged_vcf = tempfile.NamedTemporaryFile(
        prefix="unsorted", suffix=".vcf", dir=temp_dir, delete=False
    )
    concat_vcfs(split_flagged_files, unsorted_merged_vcf.name)

    unzipped_out_vcf = out_vcf.strip(".gz")
    fout = open(unzipped_out_vcf, "w")
    subprocess.check_call(["vcf-sort", unsorted_merged_vcf.name], stdout=fout)
    fout.close()

    subprocess.check_call(["bgzip", "-f", unzipped_out_vcf])
    subprocess.check_call(["tabix", "-f", "-p", "vcf", out_vcf])

    click.echo("removing temp_dir...")
    shutil.rmtree(temp_dir, ignore_errors=True)
