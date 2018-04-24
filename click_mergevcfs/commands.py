import os
import subprocess
import click
from shutil import copyfile

from click_mergevcfs.utils import get_caller, parse_header, tra2bnd

def merge_snvs(vcf_list, out_file):
    """For merging snvs and indels."""
    outdir = os.path.dirname(out_file)

    # copy input vcf to outdirs
    outdir_vcf_list = []
    for vcf in vcf_list:
        vcf_basename = os.path.basename(vcf)
        # TODO if input directory and working directory are the same, shutil.copyfile would throw a "same file" error
        copyfile(vcf, os.path.join(outdir, vcf_basename))
        outdir_vcf_list.append(os.path.join(outdir, vcf_basename))

    cmd = ["vcf-merge"]

    callers = []
    for vcf in outdir_vcf_list:
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

    fout = open(out_file, 'w')
    subprocess.check_call(cmd, stdout=fout)
    fout.close()

    parse_header(out_file, callers)

    subprocess.check_call(['bgzip', out_file])


def merge_svs(vcf_list, out_file, reference):
    """For merging svs."""
    outdir = os.path.dirname(out_file)

    # copy input vcf to outdirs
    outdir_vcf_list = []
    for vcf in vcf_list:
        vcf_basename = os.path.basename(vcf)
        copyfile(vcf, os.path.join(outdir, vcf_basename))
        outdir_vcf_list.append(os.path.join(outdir, vcf_basename))

    cmd = ["vcf-merge"]
    callers = []
    for vcf in outdir_vcf_list:
        callers.append(get_caller(vcf))
        out_vcf = vcf.split('vcf')[0] + "bnd.vcf.gz"
        tra2bnd(in_vcf=vcf, out_vcf=out_vcf, reference=reference)
        # Freshly index vcf just in case index file is older than vcf
        subprocess.check_call(['tabix', '-f', '-p', 'vcf', out_vcf])
        cmd.extend([out_vcf])

    fout = open(out_file, 'w')
    subprocess.check_call(cmd, stdout=fout)
    fout.close()
    
    # TODO parse output merged vcf header
    parse_header(out_file, callers)

    subprocess.check_call(['bgzip', out_file])

    
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
