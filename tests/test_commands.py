import os
import pytest
import shutil
import gzip

from click_mergevcfs import commands

from .utils import TEST, which

EXPECTED_SNP = "8\t29496762\t8421c1fa-505f-11e7-bd98-d60f982b8ef5\tT\tG"
EXPECTED_INDEL = "8\t117864952\ta73b5158-505d-11e7-be92-c33fa51be26e\tT\tTAA"
EXPECTED_SV = "1\t10\t.\tT\tT[1:10["

def test_run(tmpdir):
    outdir = tmpdir.strpath
    reference = TEST['reference']
    mutect_snvs = TEST['mutect_snvs']
    mutect_indels = TEST['mutect_indels']
    strelka_snvs = TEST['strelka_snvs']
    strelka_indels = TEST['strelka_indels']
    caveman_snvs = TEST['caveman_snvs']
    pindel_indels = TEST['pindel_indels']
    brass_svs = TEST['brass_svs']
    smoove_svs = TEST['smoove_svs']
    svaba_svs = TEST['svaba_svs']
    normal_bam = TEST['normal_bam']
    tumor_bam = TEST['tumor_bam']
    bedFileLoc = TEST['bedFileLoc']
    indelBed = TEST['indelBed']
    annoBedLoc = TEST['annoBedLoc']
    unmatchedVCFLoc = TEST['unmatchedVCFLoc']

    snvs_vcf = [caveman_snvs, mutect_snvs, strelka_snvs]
    indels_vcf = [pindel_indels, mutect_indels, strelka_indels]
    svs_vcf = [brass_svs, smoove_svs, svaba_svs]

    # Test snvs merge
    snvs_merged = os.path.join(outdir, "merged.snvs.vcf.gz")
    commands.merge_snvs(vcf_list=snvs_vcf, out_file=snvs_merged)
    
    with gzip.open(snvs_merged, "rb") as f:
        content = f.read()
        assert EXPECTED_SNP in content

    # Test snvs merge, check the program handles properly if the output file is not gzipped
    snvs_merged = os.path.join(outdir, "merged.snvs.vcf")
    commands.merge_snvs(vcf_list=snvs_vcf, out_file=snvs_merged)

    with open(snvs_merged, "r") as f:
        content = f.read()
        assert EXPECTED_SNP in content

    # Test indels merge
    indels_merged = os.path.join(outdir, "merged.indels.vcf.gz")
    commands.merge_snvs(vcf_list=indels_vcf, out_file=indels_merged)

    with gzip.open(indels_merged, "rb") as f:
        content = f.read()
        assert EXPECTED_INDEL in content
    
    # Test SVs merge
    svs_merged = os.path.join(outdir, "merged.svs.vcf.gz")
    commands.merge_svs(vcf_list=svs_vcf, out_file=svs_merged, reference=reference)
    with gzip.open(svs_merged, 'rb') as f:
        content = f.read()
        assert EXPECTED_SV in content

    # Test flagging
    perl_path = which("perl")
    print perl_path
    ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'click_mergevcfs'))
    flag_script = os.path.join(ROOT, "cgpFlagCaVEMan_debug.pl")
    flagConfig = os.path.join(ROOT, "flag.vcf.custom.config.ini")
    flagToVcfConfig = os.path.join(ROOT, "flag.to.vcf.custom.convert.ini")
    flagged_vcf = os.path.join(outdir, "merged.flagged.vcf.gz")
    reference_fai = reference + ".fai"

    commands.caveman_postprocess(
        perl_path=perl_path,
        flag_script=flag_script,
        in_vcf=snvs_merged,
        out_vcf=flagged_vcf,
        normal_bam=normal_bam,
        tumor_bam=tumor_bam,
        bedFileLoc=bedFileLoc,
        indelBed=indelBed,
        unmatchedVCFLoc=unmatchedVCFLoc,
        reference=reference_fai,
        flagConfig=flagConfig,
        flagToVcfConfig=flagToVcfConfig,
        annoBedLoc=annoBedLoc
        )
