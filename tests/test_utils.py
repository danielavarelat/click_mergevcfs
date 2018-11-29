"""click_mergevcfs utils tests."""

from os.path import join
import os
import gzip
import shutil
import tarfile

import pysam

from click_mergevcfs import utils
from tests.utils import TEST


def test_order_vcf_columns(tmpdir):
    """Test that mutect samples columns are always first Normal, then Tumor."""
    vcf_in_list = [TEST["mutect_snvs"], TEST["wrong_order_mutect_vcf"]]
    vcf_out_list = utils.order_vcf_columns(vcf_in_list, tmpdir.strpath)

    for vcf_out in vcf_out_list:
        with gzip.open(vcf_out, "r") as fin:
            lines = [line.decode("UTF-8") for line in fin.read().splitlines()]

        header_line_idx = [i for i, s in enumerate(lines) if s.startswith("#CHR")][0]
        header_line = lines[header_line_idx]
        header = lines[: (header_line_idx + 1)]

        for h in header:
            if h.startswith("##normal_sample="):
                normal_sample = h.split("=")[1]
            if h.startswith("##tumor_sample="):
                tumor_sample = h.split("=")[1]

        assert normal_sample == header_line.split("\t")[-2:-1][0]
        assert tumor_sample == header_line.split("\t")[-1:][0]

        vcf = pysam.VariantFile(vcf_out)
        for v in vcf:
            normal = list(v.samples)[0]
            tumor = list(v.samples)[1]
            normal_af = v.samples[normal]["AF"][0]
            tumor_af = v.samples[tumor]["AF"][0]
            assert normal_af < tumor_af


def test_decompose_multiallelic_record(tmpdir):
    """Test multiallelic variants are reported one row per allele."""
    in_vcf = TEST["mutect_indels"]
    out_vcf = join(str(tmpdir), "test_decompose_multiallelic.vcf.gz")
    utils.decompose_multiallelic_record(in_vcf=in_vcf, out_vcf=out_vcf)

    expected_record_1 = "8\t117868499\t.\tTCAC\tT"
    expected_record_2 = "8\t117868499\t.\tTCAC\tTCACC"
    with gzip.open(out_vcf, "rb") as fin:
        content = fin.read()
    assert expected_record_1 in content
    assert expected_record_2 in content

    vcf = pysam.VariantFile(out_vcf)
    for record in vcf:
        if record.pos == 117868499 and record.alts[0] == "T":
            assert record.samples[0]["AD"] == (565, 51)
            assert record.samples[0]["AF"][0] == 0.26100000739097595
            assert record.samples[0]["F1R2"] == (272, 19)
            assert record.samples[0]["F2R1"] == (284, 20)
        if record.pos == 117868499 and record.alts[0] == "TCACC":
            assert record.samples[0]["AD"] == (565, 1)
            # AF is higher than calculated
            # https://software.broadinstitute.org/gatk/documentation/article?id=11096
            assert record.samples[0]["AF"][0] == 0.04100000113248825
            assert record.samples[0]["F1R2"] == (272, 0)
            assert record.samples[0]["F2R1"] == (284, 0)


def test_rename_samples_headers(tmpdir):
    """Test the util to rename the columns samples with the callers name."""
    tmp_in_1 = join(str(tmpdir), "tmp1.in")
    tmp_in_2 = join(str(tmpdir), "tmp2.in")
    tmp_out_1 = join(str(tmpdir), "tmp1.out")
    tmp_out_2 = join(str(tmpdir), "tmp2.out")

    in_vcf_1 = TEST["strelka_snvs"]
    in_vcf_2 = TEST["germline_strelka_snvs"]
    shutil.copyfile(in_vcf_1, tmp_in_1)
    shutil.copyfile(in_vcf_2, tmp_in_2)

    utils.rename_samples_headers(tmp_in_1, tmp_out_1)
    utils.rename_samples_headers(tmp_in_2, tmp_out_2)

    expected_header_1 = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        "strelka_NORMAL\tstrelka_TUMOR\n"
    )
    expected_header_2 = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tstrelka_NORMAL\n"
    )

    with gzip.open(tmp_out_1 + ".gz", "r") as fin:
        content = fin.read()
        assert expected_header_1 in content
    with gzip.open(tmp_out_2 + ".gz", "r") as fin:
        content = fin.read()
        assert expected_header_2 in content


def test_force_symlink(tmpdir):
    """Test symlink util."""
    src = join(str(tmpdir), "src")
    dst = join(str(tmpdir), "dst")

    with open(src, "w") as f:
        f.write("Not empty.")

    utils.force_symlink(src, dst)
    assert os.path.islink(dst)


def test_force_symlink_overwrite(tmpdir):
    """Test overwrite symlink util."""
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
    """Test force link util."""
    src = join(str(tmpdir), "src")
    dst = join(str(tmpdir), "dst")

    with open(src, "w") as f:
        f.write("Not empty.")

    utils.force_link(src, dst)
    assert os.path.isfile(dst)
    assert not os.path.islink(dst)


def test_force_link_overwrite(tmpdir):
    """Test force link overwrite util."""
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
    """Test tar compression util."""
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
