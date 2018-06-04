"""click_mergevcfs commands tests."""

import os
import subprocess
import tempfile
import shutil

from click_mergevcfs.utils import get_caller, parse_header, is_gz_file, \
     decompose_multiallelic_record, add_PASSED_field, add_version


def merge_snvs(vcf_list, out_file, working_dir):
    """For merging snvs and indels."""
    working_dir_vcf_list = []
    for vcf in vcf_list:
        vcf_base_filename = os.path.basename(vcf)

        # decompose multiallelic records
        decomposed_vcf = os.path.join(
            working_dir,
            "decomposed_{}".format(vcf_base_filename)
        )
        decompose_multiallelic_record(in_vcf=vcf, out_vcf=decomposed_vcf)

        # add 'PASSED' field under INFO. ex. PASSED_caveman,PASSED_mutect
        PASSED_added_vcf = os.path.join(working_dir, vcf_base_filename)
        add_PASSED_field(in_vcf=decomposed_vcf, out_vcf=PASSED_added_vcf)

        working_dir_vcf_list.append(PASSED_added_vcf)

    cmd = ["vcf-merge", "--collapse", "none"]

    callers = []
    for vcf in working_dir_vcf_list:
        callers.append(get_caller(vcf))
        bgzip_vcf = ""
        if (not vcf.endswith('.gz')) and (not is_gz_file(vcf)):
            subprocess.check_call(['bgzip', '-f', vcf])
            bgzip_vcf = vcf + ".gz"
        else:
            bgzip_vcf = vcf
        # Freshly index vcf just in case index file is older than vcf
        subprocess.check_call(['tabix', '-f', '-p', 'vcf', bgzip_vcf])
        cmd.extend([vcf])

    cmd = list(map(str, cmd))
    temp = tempfile.NamedTemporaryFile(suffix=".tmp.vcf", dir=working_dir,
                                       delete=False)
    fout = open(temp.name, 'w')
    # Output of vcf-merge is not bgziped, regardless of the output filename
    subprocess.check_call(cmd, stdout=fout)
    fout.close()

    parse_header(vcf=temp.name, callers=callers)

    # vcf-merge may create multiple ALT alleles per record, we need to
    # break those alleles into multiple lines.
    decompose_multiallelic_record(in_vcf=temp.name, out_vcf=out_file)

    # add click_mergevcfs version to the header
    add_version(out_file)


def caveman_postprocess(perl_path, flag_script, in_vcf, out_vcf, normal_bam,
                        tumor_bam, bedFileLoc, indelBed, unmatchedVCFLoc,
                        reference, flagConfig, flagToVcfConfig, annoBedLoc):
    """Run caveman flagging on merged vcf."""
    cmd = [
        perl_path,
        flag_script,
        '-i', in_vcf,
        '-o', out_vcf,
        '-s', 'HUMAN',
        '-n', normal_bam,
        '-m', tumor_bam,
        '-b', bedFileLoc,
        '-g', indelBed,
        '-umv', unmatchedVCFLoc,
        '-ref', reference + ".fai", # Reference index (fai) from caveman help
        '-t', 'pulldown',
        '-c', flagConfig,
        '-v', flagToVcfConfig,
        '-ab', annoBedLoc,
        '--verbose'
    ]

    # Unicode to string
    cmd = list(map(str, cmd))

    subprocess.check_call(cmd)

    if out_vcf.endswith('.gz'):
        subprocess.check_call(['bgzip', out_vcf])
        shutil.move(out_vcf+'.gz', out_vcf)
