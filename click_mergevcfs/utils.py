"""click_mergevcfs utils."""

import os
import tarfile
import shutil
import gzip
import binascii
import pysam

def tra2bnd(in_vcf, out_vcf, reference):
    # TODO add SVCLASS={DEL,DUP,INV}
    tra_vcf = pysam.VariantFile(in_vcf, 'r')
    bnd_vcf = pysam.VariantFile(out_vcf, 'w', header=tra_vcf.header)
    fasta = pysam.FastaFile(reference)

    for tra_record in tra_vcf:
        bnd_record = tra_record.copy()
        bnd_record.ref = get_ref(fasta, tra_record.chrom, tra_record.pos)
        bnd_record.alts = tuple(
            [get_alt(tra_record.alts[0], bnd_record.ref, tra_record.chrom,
                     tra_record.pos)]
        )
        bnd_vcf.write(bnd_record)

    bnd_vcf.close()


def get_ref(reference, chrom, pos):
    return reference.fetch(chrom, start=int(pos)-1, end=int(pos))


def get_alt(alt, ref, chrom, pos):
    if alt == "<DEL>":
        return "{}[{}:{}[".format(ref, chrom, pos)
    if alt == "<DUP>":
        return "]{}:{}]{}".format(chrom, pos, ref)

    # if this a inversion
    if alt[0] == '[':
        return alt[:-1] + ref
    if alt[-1] == ']':
        return ref + alt[1:]
    return alt


def parse_header(vcf, callers):
    # TODO maybe instead of {}_Normal and {}_Tumor, do {}_{sample name}
    out_file = "out.tmp.vcf"
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    for c in callers:
        header += "\t{}_NORMAL\t{}_TUMOR".format(c, c)
    header += "\n"
    with open(out_file, "w") as fout:
        with open(vcf, "r") as fin:
            for line in fin:
                if line.startswith("#CHROM"):
                    fout.write(header)
                else:
                    fout.write(line)
    os.remove(vcf)
    shutil.move(out_file, vcf)


def get_caller(vcf):
    # TODO finda a better way to get caller
    callers = ["mutect", "strelka", "caveman", "pindel", "brass", "smoove",
               "svaba"]
    if is_gz_file(vcf):
        with gzip.open(vcf, 'r') as fin:
            content = fin.read()
            for c in callers:
                if c in content:
                    return c
    else:
        with open(vcf, 'r') as fin:
            content = fin.read()
            for c in callers:
                if c in content:
                    return c


def which(pgm):
    path = os.getenv('PATH')
    for p in path.split(os.path.pathsep):
        p = os.path.join(p, pgm)
        if os.path.exists(p) and os.access(p, os.X_OK):
            return p


def is_gz_file(f):
    with open(f, 'rb') as fin:
        return binascii.hexlify(fin.read(2)) == b'1f8b'


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
