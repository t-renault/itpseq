#!/usr/bin/env python3
"""Entry point for the command line tools"""

import os
import sys
import webbrowser
from pathlib import Path

import click

from . import __version__
from .config import *


class NaturalOrderGroup(click.Group):
    """
    Class to keep the commands in the order defined in the source code

    https://github.com/pallets/click/issues/513
    """

    def list_commands(self, ctx):
        return self.commands.keys()


def _print_logo():
    click.echo(
        f"""
          ▄ ▗▄▄▄▖▗▄▄▖    ▗▄▄▖▗▞▀▚▖ ▄▄▄▄             ▄▄▄▄
          ▄   █  ▐▌ ▐▌  ▐▌   ▐▛▀▀▘█   █            ▀█████▄
        ▄ █   █  ▐▛▀▘ ▄▄ ▝▀▚▖▝▚▄▄▖▀▄▄▄█ ▄▀▀▄    ▄▀▀ ▄█████
    ▀▄▄▀  █   █  ▐▌     ▗▄▄▞▘         █     ▀▄▄▀   ▀████▀
                                      ▀
                 Thank you for using itpseq!
                   version: {__version__}
    """
    )


@click.group(
    cls=NaturalOrderGroup,
    context_settings={'help_option_names': ['-h', '--help']},
    epilog='To get help on a specific command, run: itpseq COMMAND --help',
    invoke_without_command=True,
)
@click.version_option(__version__, '--version', '-v', prog_name='itpseq')
@click.pass_context
def cli(ctx):
    """itpseq - A command-line tool for the analysis of iTP-Seq data."""
    if ctx.invoked_subcommand not in {'completion', 'format'}:
        _print_logo()
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@cli.command()
def help():
    """Open the HTML documentation in the default web browser."""

    module_dir = Path(__file__).resolve().parent.parent
    docs_index = module_dir.joinpath('docs', '_build', 'html', 'index.html')

    if docs_index.exists():
        url = f'file://{docs_index}'
    else:
        print('Local documentation not found.')
        url = 'https://itpseq.readthedocs.io'
    print(f'Loading documentation from: {url}')
    webbrowser.open(url)


@cli.command(
    context_settings={'ignore_unknown_options': True}, add_help_option=False
)
@click.argument('args', nargs=-1)
def parse(args):
    """
    Parse and filter the assembled iTP-Seq FASTQ files to produce files needed for the analysis.
    """
    from . import parsing

    # Pass arguments directly to parsing.cli()
    sys.argv = ['parse'] + list(args)  # Simulate sys.argv for argparse
    parsing.main()


@cli.command()
@click.argument(
    'directory', type=click.Path(exists=True, file_okay=False, dir_okay=True)
)
@click.option(
    '--ref-labels',
    type=str,
    default='noa',
    help='Reference labels. Default is "noa".',
)
@click.option(
    '--file-pattern',
    type=str,
    default=None,
    help='Regex pattern to match the sample files in the data directory.',
)
@click.option(
    '-k',
    '--keys',
    type=str,
    default=None,
    help='Properties in the file name to use for identifying the reference.',
)
@click.option(
    '-o',
    '--output',
    type=click.Path(),
    default='itpseq_report.pdf',
    help='Output file for the report. Defaults to "itpseq_report.pdf".',
)
def report(directory, keys, ref_labels, file_pattern, output):
    """
    Generate a report for the dataset in the specified DIRECTORY.
    """
    import matplotlib

    from .core import DataSet

    matplotlib.use('Agg')

    # Create the DataSet instance
    dataset = DataSet(
        data_path=Path(directory),
        keys=tuple(keys.split(',')) if keys else keys,
        ref_labels=ref_labels,
        file_pattern=file_pattern,
    )

    # Generate the report
    click.echo(
        f'Generating report for "{directory}", please wait this may take a while...'
    )
    dataset.report(output=output)
    click.echo(f'"{output}" written')


def _files_or_dir(ctx, param, value):
    """Click callback to accept one or several files or a directory."""
    if not value:
        raise click.BadParameter(
            'At least one input file or directory must be provided.'
        )
    paths = [Path(item) for item in value]
    if len(paths) == 1 and paths[0].is_dir():
        return paths[0]
    if all(path.is_file() for path in paths):
        return paths
    raise click.BadParameter(
        'Please provide either one or more input files or a single directory containing them'
    )


@cli.command()
@click.argument(
    'files_or_dir',
    nargs=-1,
    type=click.Path(exists=True),
    callback=_files_or_dir,
)
@click.option(
    '-c', '--codons', is_flag=True, help='Group the nucleotides by codon'
)
@click.option(
    '-a', '--aa', is_flag=True, help='Display the translated peptide'
)
@click.option(
    '-r',
    '--repeat-header',
    type=int,
    default=0,
    help='Repeat header every <INTEGER> lines',
)
@click.option(
    '-o',
    '--out',
    type=click.Path(file_okay=False, dir_okay=True),
    help='Output directory',
)
@click.option(
    '-l',
    '--limit',
    type=int,
    default=0,
    help='Limit processing to <INTEGER> sequences',
)
def format(files_or_dir, codons, aa, repeat_header, limit, out):
    """
    Format the inverse-toeprint output files to display codons, amino-acids, ...

    FILES_OR_DIR: Inverse toeprint nucleotide files (ending in '.nuc.itp.txt') or directory containing them.

    \b
    Example:
    $ itpseq format nnn15_noa1.nuc.itp.txt --codons --aa --limit 5
    #                           [E] [P] [A]
                                ATG GGA CGC cccgcagtatct
                                 M   G   R  3
            ATG AGT TAC AAA GGC AAC TCG GAA caggtagcatatc
             M   S   Y   K   G   N   S   E  8
                                ATG GAA GAG gcccatgccattcc
                                 M   E   E  3
                                    ATG AAT cgaaacatgttt
                                     M   N  2
    ATG ACT ATG TTT CTT GGA CAC ACA TAA GGG aactagttaggg
     M   T   M   F   L   G   H   T   *   G  10
    """
    from .parsing import format_sequences

    format_sequences(
        files_or_dir,
        codons=codons,
        aa=aa,
        repeat_header=repeat_header,
        limit=limit,
        out=out,
    )


@cli.command()
@click.argument('shell', required=False)
def completion(shell):
    """
    Generates a completion script for the specified shell.

    To add completion for itpseq, edit your shell configuration file and add:

    \b
    bash (~/.bashrc):
      source <(itpseq completion bash)
    zsh (~/.zshenv):
      source <(itpseq completion zsh)
    fish (~/.config/fish/config.fish):
      source <(itpseq completion fish)

    Then restart your shell.
    """
    if not shell:
        # trying to detect the shell
        shell = Path(os.readlink(f'/proc/{os.getppid()}/exe')).name
        if not shell:
            raise ValueError(
                'Could not detect shell, please specify it manually.'
            )
    shell_script = os.popen(f'_ITPSEQ_COMPLETE={shell}_source itpseq').read()
    print(shell_script)


if __name__ == '__main__':
    cli()
