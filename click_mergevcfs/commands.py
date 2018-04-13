import os
import subprocess
import click

from click_mergevcfs.utils import get_caller, parse_header, tra2bnd

def merge_vcfs(vcf_list, out_file):
    """For merging snvs and indels."""
    # Freshly index vcf just in case index file is older than vcf
    cmd = ["vcf-merge"]

    callers = []
    for vcf in vcf_list:
        callers.append(get_caller(vcf))
        if not vcf.endswith('.gz'):
            subprocess.check_call(['bgzip', vcf])
        # TODO need to be indexed
        subprocess.check_call(['tabix', '-f', '-p', 'vcf', vcf])
        cmd.extend([vcf])

    fout = open(out_file, 'w')
    subprocess.check_call(cmd, stdout=fout)
    fout.close()

    parse_header(out_file, callers)


def merge_svs(vcf_list, out_file, reference):
    cmd = ['vcf-merge']

    for vcf in vcf_list:
        tra2bnd(vcf, reference)
        if not vcf.endswith('.gz'):
            subprocess.check_call(['bgzip', vcf])
        # TODO need to be indexed
        subprocess.check_call(['tabix', '-f', '-p', 'vcf', vcf])
        cmd.extend([vcf])

    fout = open(out_file, 'w')
    subprocess.check_call(cmd, stdout=fout)
    fout.close()
    

# def post_process(merged_vcf):
