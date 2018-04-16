import os
import pytest

from click_mergevcfs import commands

from .utils import TEST

EXPECTED_SNP = "8\t29496762\t8421c1fa-505f-11e7-bd98-d60f982b8ef5\tT\tG"
EXPECTED_INDEL = "8\t117864952\ta73b5158-505d-11e7-be92-c33fa51be26e\tT\tTAA"

def test_run(tmpdir):
    outdir = tmpdir.strpath
    mutect_snvs = TEST['mutect_snvs']
    mutect_indels = TEST['mutect_indels']
    strelka_snvs = TEST['strelka_snvs']
    strelka_indels = TEST['strelka_indels']
    caveman_snvs = TEST['caveman_snvs']
    pindel_indels = TEST['pindel_indels']
    # brass_sv = TEST['brass_sv']
    # svaba_sv = TEST['svaba_sv']
    # gridss_sv = TEST['gridss_sv']

    # vcfs = [
    #     mutect_snvs,
    #     mutect_indels,
    #     strelka_snvs,
    #     strelka_indels,
    #     caveman_snvs,
    #     pindel_indels,
    #     brass_sv,
    #     svaba_sv,
    #     gridss_sv,
    # ]

    snps_vcf = [caveman_snvs, mutect_snvs, strelka_snvs]
    indels_vcf = [pindel_indels, mutect_indels, strelka_indels]

    # Test snps merge
    snps_merged = os.path.join(outdir, "merged.vcf")
    commands.merge_snvs(vcf_list=snps_vcf, out_file=snps_merged)

    with open(snps_merged, "r") as f:
        content = f.read()
        assert EXPECTED_SNP in content

    # Test indels merge
    indels_merged = os.path.join(outdir, "merged.vcf")
    commands.merge_snvs(vcf_list=indels_vcf, out_file=indels_merged)

    with open(indels_merged, "r") as f:
        content = f.read()
        assert EXPECTED_INDEL in content
    
    # Test SVs merge
    # commands.merge_svs(vcf_list=svs_vcf, out_dir)
    # with open()
    # assert sum(1 for line in open(snps_merged)) == 926
    # assert sum(1 for line in open(indels_merged)) == 334