from os.path import join, dirname, abspath
import os

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
    "brass_svs": join(DATA, "svs", "brass.vcf"),
    "svaba_svs": join(DATA, "svs", "svaba.vcf"),
    "smoove_svs": join(DATA, "svs", "smoove.vcf"),
    "normal_bam": join(DATA, "bam", "normal.bam"),
    "tumor_bam": join(DATA, "bam", "tumor.bam"),
    "bedFileLoc": join(DATA, "flag_config"),
    "indelBed": join(DATA, "indel.germline.bed"),
    "annoBedLoc": join(DATA, "annotable_region"),
    "unmatchedVCFLoc": join(DATA, "pon"),
}


def which(pgm):
    path = os.getenv('PATH')
    for p in path.split(os.path.pathsep):
        p = os.path.join(p, pgm)
        if os.path.exists(p) and os.access(p, os.X_OK):
            return p
