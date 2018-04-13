"""click_mergevcfs utils tests."""

from os.path import join
import os
import tarfile

from click_mergevcfs import utils

def test_parse_header(tmpdir):
    tmp_out = join(str(tmpdir), "tmp.out")
    with open(tmp_out, 'w') as fout:
        fout.write("#CHROM\tPOS\tID")

    callers = ['caller1', 'caller2']
    utils.parse_header(tmp_out, callers)
    expected_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcaller1_NORMAL\tcaller1_TUMOR\tcaller2_NORMAL\tcaller2_TUMOR"

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
