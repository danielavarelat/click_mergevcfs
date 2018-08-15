from os.path import join, dirname, abspath

ROOT = abspath(dirname(__file__))

DATA = join(ROOT, "data")

TEST = {
    "reference": join(DATA, "reference", "reference.fasta"),
    "mutect_snvs": join(DATA, "snvs", "mutect.snvs.vcf.gz"),
    "mutect_indels": join(DATA, "indels", "mutect.indels.vcf.gz"),
    "strelka_snvs": join(DATA, "snvs", "strelka.snvs.vcf.gz"),
    "strelka_indels": join(DATA, "indels", "strelka.indels.vcf.gz"),
    "caveman_snvs": join(DATA, "snvs", "caveman.snvs.vcf.gz"),
    "pindel_indels": join(DATA, "indels", "pindel.indels.vcf.gz"),
    "normal_bam": join(DATA, "bam", "normal.bam"),
    "tumor_bam": join(DATA, "bam", "tumor.bam"),
    "bedFileLoc": join(DATA, "flag_config"),
    "indelBed": join(DATA, "indel.germline.bed"),
    "annoBedLoc": join(DATA, "annotable_region"),
    "unmatchedVCFLoc": join(DATA, "pon"),
    "wrong_order_mutect_vcf": join(DATA, "snvs", 
                                   "test_order_mutect_samples.vcf.gz"),
    "expected_vcf": join(DATA, "snvs", "expected.merged.flagged.snv.vcf.gz"),
}
