import os
import subprocess

from shutil import copyfile

from click_mergevcfs.utils import get_caller, parse_header, tra2bnd, is_gz_file, decompose_multiallelic_record, add_PASSED_field

def merge_snvs(vcf_list, out_file, working_dir):
    """
    For merging snvs and indels.
    Output file format is determined based on the out_file filename.
    """
    # copy input vcf to working_dir
    working_dir_vcf_list = []
    for vcf in vcf_list:
        vcf_base_filename = os.path.basename(vcf)

        # decompose multiallelic records
        decomposed_vcf = os.path.join(working_dir, "decomposed_{}".format(vcf_base_filename))
        decompose_multiallelic_record(in_vcf=vcf, out_vcf=decomposed_vcf)

        # add 'PASSED' field under INFO. ex. PASSED=caveman,mutect
        PASSED_added_vcf = os.path.join(working_dir, vcf_base_filename)
        add_PASSED_field(in_vcf=decomposed_vcf, out_vcf=PASSED_added_vcf)

        working_dir_vcf_list.append(PASSED_added_vcf)

    cmd = ["vcf-merge"]

    callers = []
    for vcf in working_dir_vcf_list:
        callers.append(get_caller(vcf))
        bgzip_vcf = ""
        if not vcf.endswith('.gz'):
            subprocess.check_call(['bgzip', vcf])
            bgzip_vcf = vcf + ".gz"
        else:
            bgzip_vcf = vcf
        # Freshly index vcf just in case index file is older than vcf
        subprocess.check_call(['tabix', '-f', '-p', 'vcf', bgzip_vcf])
        cmd.extend([vcf])

    cmd = list(map(str, cmd))
    fout = open(out_file, 'w')
    # Output of vcf-merge is not bgziped, regardless of the output filename
    subprocess.check_call(cmd, stdout=fout)
    fout.close()

    parse_header(out_file, callers)

    # If user specify the output file should be gziped, but the outfile is
    # not gzipped, we need to gzip the outfile
    if out_file.endswith('.gz') and (not is_gz_file(out_file)):
        corrected_filename = out_file.strip('.gz')
        os.rename(out_file, corrected_filename)
        subprocess.check_call(['bgzip', corrected_filename])


def merge_svs(vcf_list, out_file, reference, working_dir):
    """For merging svs."""
    # copy input vcf to outdirs
    working_dir_vcf_list = []
    for vcf in vcf_list:
        vcf_base_filename = os.path.basename(vcf)
        copyfile(vcf, os.path.join(working_dir, vcf_base_filename))
        working_dir_vcf_list.append(os.path.join(working_dir, vcf_base_filename))

    cmd = ["vcf-merge"]
    callers = []
    for vcf in working_dir_vcf_list:
        callers.append(get_caller(vcf))
        out_vcf = vcf.split('vcf')[0] + "bnd.vcf.gz"
        tra2bnd(in_vcf=vcf, out_vcf=out_vcf, reference=reference)
        # Freshly index vcf just in case index file is older than vcf
        subprocess.check_call(['tabix', '-f', '-p', 'vcf', out_vcf])
        cmd.extend([out_vcf])

    cmd = list(map(str, cmd))
    fout = open(out_file, 'w')
    subprocess.check_call(cmd, stdout=fout)
    fout.close()

    # TODO parse output merged vcf header
    parse_header(out_file, callers)

    # If user specify the output file should be gziped, but the outfile is not
    # gzipped, we need to gzip the outfile
    if out_file.endswith('.gz') and (not is_gz_file(out_file)):
        corrected_filename = out_file.strip('.gz')
        os.rename(out_file, corrected_filename)
        subprocess.check_call(['bgzip', corrected_filename])


def caveman_postprocess(perl_path, flag_script, in_vcf, out_vcf, normal_bam,
                        tumor_bam, bedFileLoc, indelBed, unmatchedVCFLoc,
                        reference, flagConfig, flagToVcfConfig, annoBedLoc):
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
        '-ref', reference,
        '-t', 'pulldown',
        '-c', flagConfig,
        '-v', flagToVcfConfig,
        '-ab', annoBedLoc,
        '--verbose'
    ]

    # Unicode to string
    cmd = list(map(str, cmd))

    subprocess.check_call(cmd)
