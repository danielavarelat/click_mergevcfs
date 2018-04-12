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

@click.command()
@click.option(
    "--vcf",
    multiple=True,
    help="Input vcf file"
)
@click.option(
    "--outdir",
    required=True,
    help="Path to the output directory",
    )
@click.option(
    "--no_flag",
    is_flag=True,
    default="False",
    help="Disable Caveman Postprocessing flagging"
)
@click.version_option(version=__version__)
def main(vcf, outdir, no_flag):
    """Echo message and exit."""
    commands.merge_vcfs(
        vcf_list=vcf,
        out_dir=outdir,
        no_flag=no_flag
    )

if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
