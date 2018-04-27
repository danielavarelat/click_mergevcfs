"""click_mergevcfs utils."""

import os
import tarfile
import shutil
import gzip
import binascii
import pysam

def get_caller(in_vcf):
    vcf = None
    if is_gz_file(in_vcf):
        vcf = gzip.open(in_vcf, 'rb')
    else:
        vcf = open(in_vcf, 'r')
    # Potential problem if file is too big
    content = vcf.read().lower()
    # Caveman files may contain 'pindel.germline.bed', Temporary fix
    content = content.replace('pindel.germline.bed','pindl.germline.bed')

    callers = ['mutect', 'strelka', 'mutect', 'caveman', 'pindel', 'brass']
    caller = set([c for c in callers if c in content])
    
    if len(caller) != 1:
        print "Unable to determine caller of a vcf: None or 1+ callers found."
        return False
    print "{} is produced by {}".format(in_vcf, list(caller)[0])
    return list(caller)[0]


def add_PASSED_field(in_vcf, out_vcf):
    # PASSED as string or flags?
    # see logic of merging INFO fields 
    # https: // github.com/vcftools/vcftools/blob/490848f7865abbb4b436ca09381ea7912a363fe3/src/perl/vcf-merge  # L441
    i_vcf = pysam.VariantFile(in_vcf, 'r')
    new_header = i_vcf.header.copy()
    new_header.info.add('PASSED', '.', 'String', "this variants passed which caller(s)")
    i_vcf.header.info.add('PASSED', '.', 'String', "this variants passed which caller(s)")

    o_vcf = pysam.VariantFile(out_vcf, 'w', header=new_header)

    caller = get_caller(in_vcf)

    for record in i_vcf:
        new_rec = record.copy()
        if record.filter[0] == 'PASS':
            #record.info['PASSED'] = caller
            new_rec.info['PASSED'] = caller
        o_vcf.write(new_rec)

    o_vcf.close()


def decompose_multiallelic_record(in_vcf, out_vcf):
    i_vcf = pysam.VariantFile(in_vcf, 'r')
    o_vcf = pysam.VariantFile(out_vcf, 'w', header=i_vcf.header)

    for record in i_vcf:
        number_events = len(record.alts)
        if number_events > 1:
            for i in range(0, number_events):
                new_rec = record.copy()
                new_rec.alts = tuple([record.alts[i]])
                # Mutliallic sites GT are ex. 0/1/2, which causes error later
                # Needs to change to ./.
                genotypes = list(record.samples)
                for g in genotypes:
                    new_rec.samples[g]['GT'] = (None, None)
                o_vcf.write(new_rec)
        else:
            o_vcf.write(record)

    o_vcf.close()


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
