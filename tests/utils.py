from os.path import join, dirname, abspath

ROOT = abspath(dirname(__file__))

DATA = join(ROOT, "data")

TEST = {
    "mutect_snvs": join(DATA, "snvs", "mutect.snvs.vcf.gz"),
    "mutect_indels": join(DATA, "indels", "mutect.indels.vcf.gz"),
    "strelka_snvs": join(DATA, "snvs", "strelka.snvs.vcf.gz"),
    "strelka_indels": join(DATA, "indels", "strelka.indels.vcf.gz"),
    "caveman_snvs": join(DATA, "snvs", "caveman.snvs.vcf.gz"),
    "pindel_indels": join(DATA, "indels", "pindel.indels.vcf.gz"),
    # "brass_sv": join(DATA, "sv", "brass.vcf.gz"),
    # "svaba_sv": join(DATA, "sv", "svaba.vcf.gz"),
    # "gridss_sv": join(DATA, "sv", "gridss.vcf.gz"),
}
