"""click_mergevcfs utils tests."""

from os.path import join
import os
import tarfile
import gzip
import pysam

from click_mergevcfs import utils
from .utils import TEST

# def test_tra2bnd(tmpdir):
# TODO find a variant that is both in smoove and another tra caller.
# run tra2bnd on smoove and compare
# TODO test if run tra2bnd on bnd vcf, the result is valid


def test_decompose_multiallelic_record(tmpdir):
    in_vcf = TEST['mutect_indels']
    out_vcf = join(str(tmpdir), "test_decompose_multiallelic.vcf.gz")
    utils.decompose_multiallelic_record(in_vcf=in_vcf, out_vcf=out_vcf)

    expected_record_1 = "8\t117868499\t.\tTCAC\tT"
    expected_record_2 = "8\t117868499\t.\tTCAC\tTCACC"
    with gzip.open(out_vcf, 'rb') as fin:
        content = fin.read()
    assert expected_record_1 in content
    assert expected_record_2 in content

    vcf = pysam.VariantFile(out_vcf)
    for record in vcf:
        if record.pos == 117868499 and record.alts[0] == 'T':
            assert record.samples[0]['AD'] == (565, 51)
            assert record.samples[0]['AF'][0] == 0.26100000739097595
            assert record.samples[0]['F1R2'] == (272, 19)
            assert record.samples[0]['F2R1'] == (284, 20)
        if record.pos == 117868499 and record.alts[0] == 'TCACC':
            assert record.samples[0]['AD'] == (565, 1)
            # AF is higher than calculated
            # https://software.broadinstitute.org/gatk/documentation/article?id=11096
            assert record.samples[0]['AF'][0] == 0.04100000113248825
            assert record.samples[0]['F1R2'] == (272, 0)
            assert record.samples[0]['F2R1'] == (284, 0)


def test_parse_header(tmpdir):
    tmp_out = join(str(tmpdir), "tmp.out")
    with open(tmp_out, 'w') as fout:
        fout.write("#CHROM\tPOS\tID")

    callers = ['caller1', 'caller2']
    utils.parse_header(tmp_out, callers)
    expected_header = ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                       "caller1_NORMAL\tcaller1_TUMOR\tcaller2_NORMAL\t"
                       "caller2_TUMOR")

    with open(tmp_out, 'r') as fin:
        content = fin.read()
    assert expected_header in content


def test_force_symlink(tmpdir):
    src = join(str(tmpdir), "src")
    dst = join(str(tmpdir), "dst")

    with open(src, "w") as f:
        f.write("Not empty.")

    utils.force_symlink(src, dst)
    assert os.path.islink(dst)


def test_force_symlink_overwrite(tmpdir):
    src = join(str(tmpdir), "src")
    dst = join(str(tmpdir), "dst")

    with open(src, "w") as f:
        f.write("Correct.")

    with open(dst, "w") as f:
        f.write("Wrong.")

    utils.force_symlink(src, dst)
    assert os.path.islink(dst)

    with open(dst, "r") as f:
        assert "Correct" in f.read()


def test_force_link(tmpdir):
    src = join(str(tmpdir), "src")
    dst = join(str(tmpdir), "dst")

    with open(src, "w") as f:
        f.write("Not empty.")

    utils.force_link(src, dst)
    assert os.path.isfile(dst)
    assert not os.path.islink(dst)


def test_force_link_overwrite(tmpdir):
    src = join(str(tmpdir), "src")
    dst = join(str(tmpdir), "dst")

    with open(src, "w") as f:
        f.write("Correct.")

    with open(dst, "w") as f:
        f.write("Wrong.")

    utils.force_link(src, dst)
    assert os.path.isfile(dst)
    assert not os.path.islink(dst)

    with open(dst, "r") as f:
        assert "Correct" in f.read()


def test_tar_dir(tmpdir):
    dst_dir = join(str(tmpdir), "dst_dir")
    source_dir = join(str(tmpdir), "source_dir")
    output_path = join(str(tmpdir), "source_dir.tar.gz")

    files = (
        (join(source_dir, "1"), "first file."),
        (join(source_dir, "2"), "second file."),
        (join(source_dir, "3"), "third file."),
        )

    os.makedirs(source_dir)
    os.makedirs(dst_dir)

    for i, j in files:
        with open(i, "w") as f:
            f.write(j)

    utils.tar_dir(output_path=output_path, source_dir=source_dir)
    assert tarfile.is_tarfile(output_path)

    # Remove files and dir
    for i in files:
        os.unlink(i[0])

    os.rmdir(source_dir)
    tar = tarfile.open(output_path)
    tar.extractall(path=dst_dir)
    tar.close()

    for i, j in files:
        i = i.replace(source_dir, join(dst_dir, "source_dir"))
        with open(i, "r") as f:
            assert j in f.read()
