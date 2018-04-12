import os
import subprocess
import click

def merge_vcfs(vcf_list, out_dir):
    out_file = os.path.join(out_dir, "merged.vcf")
    cmd = ["vcf-merge"]

    for vcf in vcf_list:
        cmd.extend([vcf])

    fout = open(out_file, "w")

    subprocess.check_call(cmd, stdout=fout)


# def post_process(merged_vcf):
