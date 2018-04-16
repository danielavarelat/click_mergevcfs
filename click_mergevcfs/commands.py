import os
import subprocess
import click

from click_mergevcfs.utils import get_caller, parse_header, tra2bnd

def merge_snvs(vcf_list, out_file):
    # TODO create outdir based on out_file and put new files in the outdir 
    """For merging snvs and indels."""
    outdir = os.path.dirname(out_file)

    cmd = ["vcf-merge"]

    callers = []
    for vcf in vcf_list:
        callers.append(get_caller(vcf))
        if not vcf.endswith('.gz'):
            with open(os.path.join(outdir, vcf + ".gz"), 'wb') as fout:
                subprocess.Popen(['bgzip', '-c', vcf], stdout=fout)   
        # Freshly index vcf just in case index file is older than vcf
        with open(os.path.join(outdir, vcf.split('vcf')[0] + "vcf.gz.tbi"), 'wb') as fout:
            subprocess.Popen(['tabix', '-f', '-p', 'vcf', vcf], stdout=fout)
        cmd.extend([vcf])

    fout = open(out_file, 'w')
    subprocess.check_call(cmd, stdout=fout)
    fout.close()

    parse_header(out_file, callers)


def merge_svs(vcf_list, out_file, reference):
    cmd = ['vcf-merge']

    for vcf in vcf_list:
        out_vcf = vcf.split('vcf')[0] + "bnd.vcf.gz"
        tra2bnd(in_vcf=vcf, out_vcf=out_vcf, reference=reference)
        if not vcf.endswith('.gz'):
            subprocess.check_call(['bgzip', vcf])
        # TODO need to be indexed
        subprocess.check_call(['tabix', '-f', '-p', 'vcf', vcf])
        cmd.extend([vcf])

    fout = open(out_file, 'w')
    subprocess.check_call(cmd, stdout=fout)
    fout.close()
    

# def post_process(merged_vcf):
