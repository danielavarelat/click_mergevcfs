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

from click_mergevcfs import __version__
from click_mergevcfs import commands
from click_mergevcfs import exceptions

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
    required=False,
    default="True",
    help="The input vcfs contain snvs or indels",
    )
@click.option(
    "--sv",
    required=False,
    default="False",
    help="The input vcfs contain svs",
    )
@click.option(
    "--reference",
    required=True,
    help="Genome reference file (ex. GRCH37D5)"
)
@click.option(
    "--no_flag",
    is_flag=True,
    default="False",
    help="Disable Caveman Postprocessing flagging"
)
@click.version_option(version=__version__)
def main(vcf, out, snv, sv, reference, no_flag):
    if snv and sv:
        msg = "ERROR: --snv and --sv cannot be used at the same time."
        raise exceptions.AmbiguousVariantTypeException(msg)
    if snv:
        commands.merge_snvs(vcf_list=vcf, out_file=out)
    if sv:
        commands.merge_svs(vcf_list=vcf, out_file=out, reference=reference)

if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
