"""click_mergevcfs utils."""
import os
from os.path import basename, join
import tarfile
import shutil
import gzip
import binascii
import subprocess
import tempfile

import click
import pandas as pd
from pysam import VariantFile

from click_mergevcfs import __version__


def order_vcf_columns(vcf_list, working_dir):
    """Order vcf samples if the order is not NORMAL TUMOR."""
    for index, vcf in enumerate(vcf_list):
        if get_caller(vcf) == "mutect":
            vcf_filename = basename(vcf)
            src = _order_mutect_samples(vcf)
            dst = join(working_dir, vcf_filename)
            shutil.copyfile(src, dst)
            temp = list(vcf_list)
            temp[index] = dst
            vcf_list = list(temp)
    return vcf_list


def _switch_last_two_order(list_of_string):
    return list_of_string[:-2] + list_of_string[-1:] + list_of_string[-2:-1]


def _is_normal_tumor(column_names, normal_sample, tumor_sample):
    return (
        normal_sample == column_names[-2:-1][0] and tumor_sample == column_names[-1:][0]
    )


def _order_mutect_samples(in_vcf_gz):
    """Order the mutect output vcf sample columns as Normal Tumor."""
    with gzip.open(in_vcf_gz, "r") as f:
        lines = [l.decode("UTF-8") for l in f.read().splitlines()]

    # find the line number of the header
    header_line_idx = [i for i, s in enumerate(lines) if s.startswith("#CHR")][0]

    header = lines[: (header_line_idx + 1)]
    old_df = pd.read_csv(in_vcf_gz, sep="\t", header=header_line_idx)

    column_names = list(old_df)
    for h in header:
        if h.startswith("##normal_sample="):
            normal_sample = h.split("=")[1]
        if h.startswith("##tumor_sample="):
            tumor_sample = h.split("=")[1]

    if not _is_normal_tumor(column_names, normal_sample, tumor_sample):
        click.echo(
            "tumor_sample {} and normal_sample {} are not in N T "
            "order in {}".format(tumor_sample, normal_sample, in_vcf_gz)
        )
        # switch the last two columns
        new_column_order = _switch_last_two_order(column_names)
        new_df = old_df[new_column_order]

        # change the order in the header
        header[-1] = "\t".join(_switch_last_two_order(header[-1:][0].split("\t")))

        temp_file = tempfile.mkstemp()[1]
        with open(temp_file, "w") as f:
            for l in header:
                f.write(l + "\n")
            new_df.to_csv(f, sep="\t", index=False, header=False)

        subprocess.check_call(["bgzip", temp_file])

        return temp_file + ".gz"

    return in_vcf_gz


def get_caller(in_vcf):
    """Determine which caller produced a given vcf file."""
    vcf = None
    if is_gz_file(in_vcf):
        vcf = gzip.open(in_vcf, "rb")
    else:
        vcf = open(in_vcf, "r")
    # Potential problem if file is too big
    content = vcf.read().lower()
    # Caveman files may contain 'pindel.germline.bed', Temporary fix
    content = content.replace(b"pindel.germline.bed", b"pindl.germline.bed")

    callers = ["mutect", "strelka", "caveman", "pindel", "freebayes"]
    caller = {c for c in callers if c.encode("utf-8") in content}

    if len(caller) != 1:
        sample_name = os.path.basename(in_vcf).split(".vcf")[0]
        click.echo(
            "Unable to determine caller of {}: None or 1+ callers found. "
            "Using sample name {}."
        ).format(in_vcf, sample_name)
        return sample_name
    click.echo("{} is produced by {}".format(in_vcf, list(caller)[0]))
    return list(caller)[0]


def add_PASSED_field(in_vcf, out_vcf):
    """
    Add PASSED_{caller} fields.

    Add flags (e.g. PASSED_caveman) under INFO for PASS variant in aim of reduce
    ambiguity of confident variants in the merged vcf.
    """
    # see logic of merging INFO fields
    # https://github.com/vcftools/vcftools/blob/490848f7865abbb4b436ca09381ea7912a363fe3/src/perl/vcf-merge
    caller = get_caller(in_vcf)

    i_vcf = VariantFile(in_vcf, "rb")
    new_header = i_vcf.header.copy()
    try:
        new_header.info.add(
            "PASSED_{}".format(caller),
            ".",
            "Flag",
            "this variants passed which caller(s)",
        )
        i_vcf.header.info.add(
            "PASSED_{}".format(caller),
            ".",
            "Flag",
            "this variants passed which caller(s)",
        )
    except ValueError:
        pass

    raw_out = out_vcf.strip(".gz")
    o_vcf = VariantFile(raw_out, "w", header=new_header)

    for record in i_vcf:
        new_rec = record.copy()
        filters = list(record.filter)
        if filters and filters[0] == "PASS":
            new_rec.info["PASSED_{}".format(caller)] = 1
        o_vcf.write(new_rec)

    o_vcf.close()

    subprocess.check_call(["bgzip", "-f", raw_out])


def decompose_multiallelic_record(in_vcf, out_vcf):
    """Break records with multiple ALT alleles into multiple records."""
    i_vcf = VariantFile(in_vcf, "r")
    raw_out = out_vcf.strip(".gz")
    o_vcf = VariantFile(raw_out, "w", header=i_vcf.header)

    for record in i_vcf:
        # Only mutect put multiple ALTs in one record
        number_events = len(record.alts)
        # Temporary fix due to segfault
        # see https://github.com/leukgen/click_mergevcfs/issues/2
        if number_events >= 8:
            continue
        elif number_events > 1:
            click.echo("file={},pos={}".format(in_vcf, record.pos))
            for i in range(0, number_events):
                new_rec = record.copy()
                new_rec.alts = tuple([record.alts[i]])
                # Multiallic sites GT are ex. 0/1/2, which causes error later
                # Needs to change to ./.
                genotypes = list(record.samples)
                for g in genotypes:
                    # Overwrite GT
                    new_rec.samples[g]["GT"] = (None, None)
                    # Use none_if_tuple_out_of_idx because
                    # record.samples[g]['AD'] would sometimes return
                    # a tuple of (None,)
                    if "AD" in list(record.samples[g]):
                        new_rec.samples[g]["AD"] = (
                            record.samples[g]["AD"][0],
                            none_if_tuple_out_of_idx(
                                t=record.samples[g]["AD"], index=i + 1
                            ),
                        )
                    if "AF" in list(record.samples[g]):
                        new_rec.samples[g]["AF"] = none_if_tuple_out_of_idx(
                            t=record.samples[g]["AF"], index=i
                        )
                    if "F1R2" in list(record.samples[g]):
                        new_rec.samples[g]["F1R2"] = (
                            record.samples[g]["F1R2"][0],
                            none_if_tuple_out_of_idx(
                                t=record.samples[g]["F1R2"], index=i + 1
                            ),
                        )
                    if "F2R1" in list(record.samples[g]):
                        new_rec.samples[g]["F2R1"] = (
                            record.samples[g]["F2R1"][0],
                            none_if_tuple_out_of_idx(
                                t=record.samples[g]["F2R1"], index=i + 1
                            ),
                        )
                o_vcf.write(new_rec)
        else:
            o_vcf.write(record)

    o_vcf.close()
    subprocess.check_call(["bgzip", "-f", raw_out])


def get_ref(reference, chrom, pos):
    """Return the reference base at a given genomic position."""
    return reference.fetch(chrom, start=int(pos) - 1, end=int(pos))


def rename_samples_headers(in_vcf, out_vcf):
    """Replace hard-to-read and ambiguious header with clear header."""
    out_vcf = out_vcf.strip(".gz")
    caller = get_caller(in_vcf)
    vcf = VariantFile(in_vcf, "rb")
    samples_names = list(vcf.header.samples)
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

    if len(samples_names) == 2:
        header += "\t{}_NORMAL\t{}_TUMOR\n".format(caller, caller)
    elif len(samples_names) == 1:
        header += "\t{}_NORMAL\n".format(caller)

    with open(out_vcf, "w") as fout:
        with gzip.open(in_vcf, "r") as fin:
            for line in fin:
                if line.startswith("#CHROM"):
                    fout.write(header)
                else:
                    fout.write(line)

    subprocess.check_call(["bgzip", "-f", out_vcf])
    os.remove(in_vcf)


def which(pgm):
    """Python equivalent of linux `which`."""
    path = os.getenv("PATH")
    for p in path.split(os.path.pathsep):
        p = os.path.join(p, pgm)
        if os.path.exists(p) and os.access(p, os.X_OK):
            return p


def add_version(in_vcf):
    """Add click_mergevcfs version in the output vcf header."""
    temp = tempfile.NamedTemporaryFile(suffix=".tmp.vcf", delete=False)
    if in_vcf.endswith("gz"):
        with gzip.open(in_vcf, "rb") as fin:
            lines = fin.readlines()
    else:
        with open(in_vcf, "r") as fin:
            lines = fin.readlines()
    lines.insert(1, "##click_mergevcfs={}\n".format(__version__))
    with open(temp.name, "w") as fout:
        for l in lines:
            fout.write(l.decode("utf-8") if isinstance(l, bytes) else l)
    subprocess.check_call(["bgzip", temp.name])
    shutil.move(temp.name + ".gz", in_vcf)


def is_gz_file(f):
    """Return true if a given file is gziped."""
    with open(f, "rb") as fin:
        return binascii.hexlify(fin.read(2)) == b"1f8b"


def none_if_tuple_out_of_idx(t, index):
    """Return none if t[i] is out of index, else return t[i]."""
    return t[index] if len(t) > index else None


def force_link(src, dst):
    """Force a link between src and dst."""
    try:
        os.unlink(dst)
        os.link(src, dst)
    except OSError:
        os.link(src, dst)


def force_symlink(src, dst):
    """Force a symlink between src and dst."""
    try:
        os.unlink(dst)
        os.symlink(src, dst)
    except OSError:
        os.symlink(src, dst)


def tar_dir(output_path, source_dir):
    """Compress a `source_dir` in `output_path`."""
    with tarfile.open(output_path, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
