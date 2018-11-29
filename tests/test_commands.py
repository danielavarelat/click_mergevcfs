"""Module to test click_mergevcfs commands."""
import os
from os.path import abspath, dirname, join
import gzip

from click_mergevcfs import commands
from click_mergevcfs import utils
from click_mergevcfs import __version__
from .utils import get_variant_lines, TEST


# Expected variants in test data
EXPECTED_SOMATIC_SNP = "8\t29496762\t8421c1fa-505f-11e7-bd98-d60f982b8ef5\tT\tG"
EXPECTED_SOMATIC_INDEL = "8\t117864952\ta73b5158-505d-11e7-be92-c33fa51be26e\tT\tTAA"

EXPECTED_GERMLINE_SNP_1 = "2\t796\t.\tT\tA"
EXPECTED_GERMLINE_SNP_2 = "2\t123881\t.\tT\tC"
EXPECTED_GERMLINE_INDEL_1 = "2\t123513\t.\tT\tTG"
EXPECTED_GERMLINE_INDEL_2 = "2\t123601\t.\tTG\tT"

# Test data
REFERENCE = TEST["reference"]
NORMAL_BAM = TEST["normal_bam"]
TUMOR_BAM = TEST["tumor_bam"]
BEDFILELOC = TEST["bedFileLoc"]
INDELBED = TEST["indelBed"]
ANNOBEDLOC = TEST["annoBedLoc"]
UNMATCHEDVCFLOC = TEST["unmatchedVCFLoc"]

CAVEMAN_SNVS = TEST["caveman_snvs"]
PINDEL_INDELS = TEST["pindel_indels"]
MUTECT_SNVS = TEST["mutect_snvs"]
MUTECT_INDELS = TEST["mutect_indels"]
STRELKA_SNVS = TEST["strelka_snvs"]
STRELKA_INDELS = TEST["strelka_indels"]
GERMLINE_STRELKA_SNVS = TEST["germline_strelka_snvs"]
GERMLINE_STRELKA_INDELS = TEST["germline_strelka_indels"]
GERMLINE_FREEBAYES_SNVS = TEST["germline_freebayes_snvs"]
GERMLINE_FREEBAYES_INDELS = TEST["germline_freebayes_indels"]


def run_merge_command(
    vcf_list,
    vcf_out,
    working_dir,
    expected_output,
    expected_count,
    not_expected_output=None,
):
    """Test the vcf snvs merge of the triple caller."""
    os.makedirs(working_dir)
    commands.merge_snvs(vcf_list=vcf_list, out_file=vcf_out, working_dir=working_dir)

    assert utils.is_gz_file(vcf_out)
    with gzip.open(vcf_out, "rb") as f:
        content = f.read()
        assert get_variant_lines(content) == expected_count
        for expected_content in expected_output:
            assert expected_content in content
        for not_expected_content in not_expected_output or []:
            assert not_expected_content not in content

    # Test if the final output only contains one ALT per record
    with gzip.open(vcf_out, "rb") as f:
        for record in f:
            if not record.startswith("#"):
                assert not "," in record.split("\t")[4]


def test_run_snvs(tmpdir):
    """Test the vcf snvs merge of the triple caller."""
    outdir = tmpdir.strpath
    vcf_list = [CAVEMAN_SNVS, MUTECT_SNVS, STRELKA_SNVS]
    vcf_out = join(outdir, "merged.snvs.vcf.gz")
    working_dir = join(outdir, "workingdir", "snv")
    expected_output = [
        "##click_mergevcfs={}\n".format(__version__),
        "caveman_NORMAL",
        "caveman_TUMOR",
        "strelka_NORMAL",
        "strelka_TUMOR",
        "mutect_NORMAL",
        "mutect_TUMOR",
        EXPECTED_SOMATIC_SNP,
        "PASSED_strelka;",
        "PASSED_mutect;",
    ]
    expected_count = 707
    run_merge_command(
        vcf_list=vcf_list,
        vcf_out=vcf_out,
        working_dir=working_dir,
        expected_output=expected_output,
        expected_count=expected_count,
    )


def test_run_indels(tmpdir):
    """Test the vcf indels merge of the triple caller."""
    outdir = tmpdir.strpath
    vcf_list = [PINDEL_INDELS, MUTECT_INDELS, STRELKA_INDELS]
    vcf_out = join(outdir, "merged.indels.vcf.gz")
    working_dir = join(outdir, "workingdir", "indels")
    expected_output = [
        "##click_mergevcfs={}\n".format(__version__),
        "pindel_NORMAL",
        "pindel_TUMOR",
        "strelka_NORMAL",
        "strelka_TUMOR",
        "mutect_NORMAL",
        "mutect_TUMOR",
        EXPECTED_SOMATIC_INDEL,
        "PASSED_pindel;",
        "PASSED_strelka;",
        "PASSED_mutect;",
    ]
    expected_count = 118
    run_merge_command(
        vcf_list=vcf_list,
        vcf_out=vcf_out,
        working_dir=working_dir,
        expected_output=expected_output,
        expected_count=expected_count,
    )


def test_run_germline(tmpdir):
    """Test the vcf merge of the germline pipeline."""
    outdir = tmpdir.strpath
    vcf_list = [GERMLINE_FREEBAYES_SNVS, GERMLINE_STRELKA_SNVS]
    vcf_out = join(outdir, "merged.germline.snvs.vcf.gz")
    working_dir = join(outdir, "workingdir", "snvs")
    expected_output = [
        "##click_mergevcfs={}\n".format(__version__),
        "freebayes_NORMAL",
        "strelka_NORMAL",
        EXPECTED_GERMLINE_SNP_1,
        EXPECTED_GERMLINE_SNP_2,
        "PASSED_strelka;",
    ]
    not_expected_output = ["freebayes_TUMOR", "strelka_TUMOR"]
    expected_count = 16
    run_merge_command(
        vcf_list=vcf_list,
        vcf_out=vcf_out,
        working_dir=working_dir,
        expected_output=expected_output,
        expected_count=expected_count,
        not_expected_output=not_expected_output,
    )

    outdir = tmpdir.strpath
    vcf_list = [GERMLINE_FREEBAYES_INDELS, GERMLINE_STRELKA_INDELS]
    vcf_out = join(outdir, "merged.germline.indels.vcf.gz")
    working_dir = join(outdir, "workingdir", "indels")
    expected_output = [
        "##click_mergevcfs={}\n".format(__version__),
        "freebayes_NORMAL",
        "strelka_NORMAL",
        EXPECTED_GERMLINE_INDEL_1,
        EXPECTED_GERMLINE_INDEL_2,
        "PASSED_strelka;",
    ]
    not_expected_output = ["freebayes_TUMOR", "strelka_TUMOR"]
    expected_count = 2
    run_merge_command(
        vcf_list=vcf_list,
        vcf_out=vcf_out,
        working_dir=working_dir,
        expected_output=expected_output,
        expected_count=expected_count,
        not_expected_output=not_expected_output,
    )


def test_run_flagging(tmpdir):
    """Test snvs flagging."""
    outdir = tmpdir.strpath
    vcf_list = [CAVEMAN_SNVS, MUTECT_SNVS, STRELKA_SNVS]
    flagging_temp_dir = join(outdir, "workingdir", "flagging")
    os.makedirs(flagging_temp_dir)
    vcf_out = join(outdir, "merged.snvs.vcf.gz")
    commands.merge_snvs(
        vcf_list=vcf_list, out_file=vcf_out, working_dir=flagging_temp_dir
    )

    # Test flagging
    perl_path = utils.which("perl")
    utils_dir = join(abspath(dirname(__file__)), os.pardir, "utils")
    flag_script = join(utils_dir, "cgpFlagCaVEMan_custom.pl")
    flag_config = join(utils_dir, "flag.vcf.custom.config.ini")
    flag_to_vcf_config = join(utils_dir, "flag.to.vcf.custom.convert.ini")
    flagged_vcf = join(outdir, "merged.flagged.snv.vcf.gz")
    expected_flaggd_vcf = TEST["expected_vcf"]

    commands.caveman_postprocess(
        perl_path=perl_path,
        flag_script=flag_script,
        in_vcf=vcf_out,
        out_vcf=flagged_vcf,
        bin_size=300,  # small bin_size for testing correctness of parallizing
        working_dir=flagging_temp_dir,
        normal_bam=NORMAL_BAM,
        tumor_bam=TUMOR_BAM,
        bedFileLoc=BEDFILELOC,
        indelBed=INDELBED,
        unmatchedVCFLoc=UNMATCHEDVCFLOC,
        reference=REFERENCE,
        flagConfig=flag_config,
        flagToVcfConfig=flag_to_vcf_config,
        annoBedLoc=ANNOBEDLOC,
        sequencing_method="TGD",
    )

    # Test if the flagged vcf is gziped
    assert utils.is_gz_file(flagged_vcf)

    # Test to see if split is correct
    with gzip.open(expected_flaggd_vcf) as f:
        lines = f.readlines()
    variants = [l for l in lines if not l.startswith("#")]
    expected_num_variants = len(variants)

    with gzip.open(flagged_vcf) as f:
        lines = f.readlines()
    variants = [l for l in lines if not l.startswith("#")]
    obs_num_variants = len(variants)

    assert expected_num_variants == obs_num_variants
