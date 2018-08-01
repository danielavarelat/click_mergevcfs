"""click_mergevcfs utils."""

import os
import tarfile
import shutil
import gzip
import binascii
import subprocess
import tempfile
import pysam
import pandas as pd

from click_mergevcfs import __version__

def order_mutect_samples(in_vcf_gz):
    """Order the mutect output vcf sample columns as Normal Tumor."""
    def switch_last_two_order(list_of_string):
        return list_of_string[:-2] + list_of_string[-1:] + list_of_string[-2:-1]

    def normal_tumor(column_names):
        if 'N' in column_names[-2:-1][0] and 'T' in column_names[-1:][0]:
            return True
        return False

    with gzip.open(in_vcf_gz, 'r') as f:
        lines = [l.decode('UTF-8') for l in f.read().splitlines()]

    # find the line number of the header
    n = 0
    for l in lines:
        if l.startswith("#CHROM"):
            break
        n += 1

    header = lines[:(n+1)]
    old_df = pd.read_csv(in_vcf_gz, sep='\t', header=n)

    column_names = list(old_df)
    if not normal_tumor(column_names):
        # switch the last two columns
        new_column_order = switch_last_two_order(column_names)
        new_df = old_df[new_column_order]

        # change the order in the header
        header[-1] = '\t'.join(
            switch_last_two_order(header[-1:][0].split('\t')))

        temp_file = tempfile.mkstemp()[1]
        with open(temp_file, 'w') as f:
            for l in header:
                f.write(l+'\n')
            new_df.to_csv(f, sep='\t', index=False, header=False)

        subprocess.check_call(['bgzip', temp_file])

        return temp_file + '.gz'

    return in_vcf_gz


def get_caller(in_vcf):
    """Determine which caller produced a given vcf file."""
    vcf = None
    if is_gz_file(in_vcf):
        vcf = gzip.open(in_vcf, 'rb')
    else:
        vcf = open(in_vcf, 'r')
    # Potential problem if file is too big
    content = vcf.read().lower()
    # Caveman files may contain 'pindel.germline.bed', Temporary fix
    content = content.replace(b'pindel.germline.bed', b'pindl.germline.bed')

    callers = ['mutect', 'strelka', 'caveman', 'pindel']
    caller = set([c for c in callers if c.encode('utf-8') in content])

    if len(caller) != 1:
        sample_name = os.path.basename(in_vcf).split(".vcf")[0]
        print ("Unable to determine caller of {}: None or 1+ callers found. "
               "Using sample name {}.").format(in_vcf, sample_name)
        return sample_name
    print("{} is produced by {}".format(in_vcf, list(caller)[0]))
    return list(caller)[0]


def add_PASSED_field(in_vcf, out_vcf):
    """
    Add PASSED_{caller} fields.

    Add flags (e.g. PASSED_caveman) under INFO for PASS variant in aim of reduce
    ambiguity of confident variants in the merged vcf.
    """
    # see logic of merging INFO fields
    # pylint: disable=C0321
    # https://github.com/vcftools/vcftools/blob/490848f7865abbb4b436ca09381ea7912a363fe3/src/perl/vcf-merge
    caller = get_caller(in_vcf)

    i_vcf = pysam.VariantFile(in_vcf, 'rb')
    new_header = i_vcf.header.copy()
    try:
        new_header.info.add('PASSED_{}'.format(caller), '.', 'Flag',
                            "this variants passed which caller(s)")
        i_vcf.header.info.add('PASSED_{}'.format(caller), '.', 'Flag',
                              "this variants passed which caller(s)")
    except ValueError:
        pass

    raw_out = out_vcf.strip('.gz')
    o_vcf = pysam.VariantFile(raw_out, 'w', header=new_header)

    for record in i_vcf:
        new_rec = record.copy()
        if list(record.filter)[0] == 'PASS':
            new_rec.info['PASSED_{}'.format(caller)] = 1
        o_vcf.write(new_rec)

    o_vcf.close()

    subprocess.check_call(['bgzip', '-f', raw_out])


def decompose_multiallelic_record(in_vcf, out_vcf):
    """Break records with multiple ALT alleles into multiple records."""
    i_vcf = pysam.VariantFile(in_vcf, 'r')
    raw_out = out_vcf.strip('.gz')
    o_vcf = pysam.VariantFile(raw_out, 'w', header=i_vcf.header)

    for record in i_vcf:
        number_events = len(record.alts)
        # Only mutect put multiple ALTS in one record
        if number_events > 1:
            print("file={},pos={}".format(in_vcf, record.pos))
            for i in range(0, number_events):
                # Temporary fix due to segfault
                # see https://github.com/leukgen/click_mergevcfs/issues/2
                if len(record.samples[0]['AF']) >= 8:
                    continue
                new_rec = record.copy()
                new_rec.alts = tuple([record.alts[i]])
                # Mutliallic sites GT are ex. 0/1/2, which causes error later
                # Needs to change to ./.
                genotypes = list(record.samples)
                for g in genotypes:
                    # Overwrite GT
                    new_rec.samples[g]['GT'] = (None, None)
                    # Use none_if_tuple_out_of_idx because
                    # record.samples[g]['AD'] would sometimes return
                    # a tuple of (None,)
                    if 'AD' in list(record.samples[g]):
                        new_rec.samples[g]['AD'] = (
                            record.samples[g]['AD'][0],
                            none_if_tuple_out_of_idx(
                                t=record.samples[g]['AD'],
                                index=i+1
                            )
                        )
                    if 'AF' in list(record.samples[g]):
                        new_rec.samples[g]['AF'] = none_if_tuple_out_of_idx(
                            t=record.samples[g]['AF'],
                            index=i
                        )
                    if 'F1R2' in list(record.samples[g]):
                        new_rec.samples[g]['F1R2'] = (
                            record.samples[g]['F1R2'][0],
                            none_if_tuple_out_of_idx(
                                t=record.samples[g]['F1R2'],
                                index=i+1
                            )
                        )
                    if 'F2R1' in list(record.samples[g]):
                        new_rec.samples[g]['F2R1'] = (
                            record.samples[g]['F2R1'][0],
                            none_if_tuple_out_of_idx(
                                t=record.samples[g]['F2R1'],
                                index=i+1
                            )
                        )
                o_vcf.write(new_rec)
        else:
            o_vcf.write(record)

    o_vcf.close()

    subprocess.check_call(['bgzip', '-f', raw_out])


def get_ref(reference, chrom, pos):
    """Return the reference base at a given genomic position."""
    return reference.fetch(chrom, start=int(pos)-1, end=int(pos))


def parse_header(vcf, callers):
    """Replace hard-to-read and ambiguious header with clear header."""
    # TODO maybe instead of {}_Normal and {}_Tumor, do {}_{sample name}
    temp = tempfile.NamedTemporaryFile()
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    for c in callers:
        header += "\t{}_NORMAL\t{}_TUMOR".format(c, c)
    header += "\n"
    with open(temp.name, "w") as fout:
        with open(vcf, "r") as fin:
            for line in fin:
                if line.startswith("#CHROM"):
                    fout.write(header)
                else:
                    fout.write(line)
    os.remove(vcf)
    shutil.copyfile(temp.name, vcf)

def which(pgm):
    """Python equivlent of linux `which`."""
    path = os.getenv('PATH')
    for p in path.split(os.path.pathsep):
        p = os.path.join(p, pgm)
        if os.path.exists(p) and os.access(p, os.X_OK):
            return p

def add_version(in_vcf):
    """Add click_mergevcfs version in the output vcf header."""
    temp = tempfile.NamedTemporaryFile(suffix=".tmp.vcf", delete=False)
    if in_vcf.endswith('gz'):
        with gzip.open(in_vcf, 'rb') as fin:
            lines = fin.readlines()
    else:
        with open(in_vcf, 'r') as fin:
            lines = fin.readlines()
    lines.insert(1, "##click_mergevcfs={}\n".format(__version__))
    with open(temp.name, 'w') as fout:
        for l in lines:
            fout.write(l.decode("utf-8") if isinstance(l, bytes) else l)
    subprocess.check_call(["bgzip", temp.name])
    shutil.move(temp.name+".gz", in_vcf)

def is_gz_file(f):
    """Return true if a given file is gziped."""
    with open(f, 'rb') as fin:
        return binascii.hexlify(fin.read(2)) == b'1f8b'

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
