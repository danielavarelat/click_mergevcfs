"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?
You might be tempted to import things from __main__ later, but that will
cause problems, the code will get executed twice:

    - When you run `python -m click_mergevcfs` python will execute
      `__main__.py` as a script. That means there won't be any
      `click_mergevcfs.__main__` in `sys.modules`.

    - When you import __main__ it will get executed again (as a module) because
      there's no `click_mergevcfs.__main__` in `sys.modules`.

Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""

import click
from os.path import join, dirname, abspath

from click_mergevcfs import __version__
from click_mergevcfs import commands
from click_mergevcfs import exceptions
from click_mergevcfs import utils

@click.command()
@click.option(
    "--vcf",
    multiple=True,
    help="Input vcf file"
)
@click.option(
    "--out",
    required=True,
    help="Path to the output file",
    )
@click.option(
    "--snv",
    is_flag=True,
    help="The input vcfs contain snvs or indels",
    )
@click.option(
    "--sv",
    is_flag=True,
    help="The input vcfs contain svs",
    )
@click.option(
    "--reference",
    help="Genome reference file (ex. GRCH37D5)"
)
@click.option(
    "--no_flag",
    is_flag=True,
    default="False",
    help="Disable Caveman Postprocessing flagging"
)
@click.version_option(version=__version__)
def main(vcf, outdir, snv, sv, reference, no_flag, **kwargs):
    if snv and sv:
        msg = "ERROR: --snv and --sv cannot be used at the same time."
        raise exceptions.AmbiguousVariantTypeException(msg)
    if snv:
        # TODO better file name with sample name
        merged_vcf = join(outdir, "merged.snv.vcf.gz")
        commands.merge_snvs(vcf_list=vcf, out_file=merged_vcf)
    if sv:
        merged_vcf = join(outdir, "merged.sv.vcf.gz")
        commands.merge_svs(vcf_list=vcf, out_file=merged_vcf, reference=reference)
    else:
        msg = "ERROR: no variant type is specified in the options. Use either --snv or --sv."
        raise exceptions.AmbiguousVariantTypeException(msg)

    if not no_flag:
        # TODO check parameter

        perl_path = utils.which('perl')
        ROOT = abspath(dirname(__file__))
        flag_script = join(ROOT, "cgpFlagCaVEMan_debug.pl")
        flagConfig = join(ROOT, "flag.vcf.custom.config.ini")
        flagToVcfConfig = join(ROOT, "flag.to.vcf.custom.convert.ini")
        flagged_vcf = join(outdir, "merged.flagged.vcf.gz")
        
        commands.caveman_postprocess(
            perl_path=perl_path,
            flag_script=flag_script,
            in_vcf=merged_vcf,
            out_vcf=flagged_vcf,
            normal_bam=kwargs['normal_bam'],
            tumor_bam=kwargs['tumor_bam'],
            bedFileLoc=kwargs['bedFileLoc'],
            indelBed=kwargs['indelBed'],
            unmatchedVCFLoc=kwargs['unmatchedVCFLoc'],
            reference=reference,
            flagConfig=flagConfig,
            flagToVcfConfig=flagToVcfConfig,
            annoBedLoc=kwargs['annoBedLoc']
            )

        
if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
