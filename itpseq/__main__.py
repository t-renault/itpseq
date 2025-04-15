#!/usr/bin/env python3
"""Entry point for the command line tools"""

import os
import sys
import webbrowser
from pathlib import Path

import rich
import rich_click as click
from rich.console import Console

from . import __version__

click.rich_click.MAX_WIDTH = 80


def arg_range(ctx, param, value):
    try:
        if '-' in value:
            start, end = value.split('-')
            start = int(start) if start else 0
            end = int(end) if end else 0
        else:
            start = end = int(value)
        return (start, end)
    except ValueError:
        raise click.BadParameter(
            'Range must be in format "min-max" (1-10), "min-" (1-), "-max" (-10), or "value" (10).'
        )


# disable verbose output of IntRange
class CustomIntRange(click.IntRange):
    def get_metavar(self, param):
        return f'{self.min}-{self.max}'

    def _describe_range(self):
        return ''

    def __repr__(self):
        return f'<{self.min}-{self.max}>'


class NaturalOrderGroup(click.RichGroup):
    """
    Class to keep the commands in the order defined in the source code

    https://github.com/pallets/click/issues/513
    """

    def list_commands(self, ctx):
        return self.commands.keys()


def _print_logo():
    VERSION = f'version: {__version__}'.center(54)
    console = Console(style='')
    console.print(
        f"""
           ▄ ▗▄▄▄▖▗▄▄▖    ▗▄▄▖▗▞▀▚▖ ▄▄▄▄            ^▄▄▄▄$
           ▄   █  ▐▌ ▐▌  ▐▌   ▐▛▀▀▘█   █           ^▀█████▄$
    £    ▄$█   █  ▐▛▀▘£▄▄$▝▀▚▖▝▚▄▄▖▀▄▄▄█£▄▀▀▄    ▄▀▀^▄█████$
    £▀▄▄▀ $█   █  ▐▌     ▗▄▄▞▘         █£    ▀▄▄▀  ^▀████▀$
                                       ▀
                  Thank you for using itpseq!
                 https://itpseq.readthedocs.io
     {                   VERSION                           }
    """.replace(
            '^', ' [grey70]'
        )
        .replace('£', ' [grey50]')
        .replace('$', ' [default]')
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


@cli.command(cls=click.RichCommand)
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
    cls=click.RichCommand,
    context_settings={'show_default': True},
    # context_settings={'ignore_unknown_options': True}, add_help_option=False
)
@click.argument(
    'files', nargs=-1, type=click.Path(exists=True, dir_okay=False)
)
@click.option(
    '--outdir',
    '-o',
    type=click.Path(dir_okay=True, file_okay=False, writable=True),
    default='.',
    help='Directory for the output files.',
)
@click.option(
    '-a1',
    '--left-adaptor',
    'a1',  # use a1 as variable name
    type=click.STRING,
    default='GTATAAGGAGGAAAAAAT',
    help='Sequence (5ʹ→3ʹ) of the left adaptor',
)
@click.option(
    '-a2',
    '--right-adaptor',
    'a2',  # use a2 as variable name
    type=click.STRING,
    default='GGTATCTCGGTGTGACTG',
    help='Sequence (5ʹ→3ʹ) of the right adaptor',
)
@click.option(
    '-s',
    '--peptide-size',
    type=click.UNPROCESSED,
    callback=arg_range,
    default='1-10',
    help='Allowed length range of peptide considered ("min-max" or "min-" or "-max" or "exact")',
)
@click.option(
    '-q',
    '--quality',
    type=CustomIntRange(min=0, max=60),
    default=30,
    help='Threshold for the PHRED quality cutoff',
)
@click.option(
    '-l',
    '--limit',
    type=click.INT,
    default=None,
    help='Max number of sequences to process (useful for quick tests)',
)
def parse(
    **kwargs,
):  # left_adaptor, right_adaptor, peptide_size, quality, limit):
    """
    Parse and filter the assembled iTP-Seq FASTQ files to produce files needed for the analysis.
    """
    import datetime

    import rich.pretty

    from .parsing import parse_all

    min_seq_len, max_seq_len = kwargs.pop('peptide_size')
    # convert codon to nucleotide positions
    kwargs['min_seq_len'] = min_seq_len * 3
    kwargs['max_seq_len'] = max_seq_len * 3 or None

    # display parameters passed to parse_all
    rich.pretty.pprint(kwargs)

    click.echo(datetime.datetime.now().strftime('Start time: %Y-%m-%d %H:%M'))
    parse_all(
        **kwargs,
        save=True,
    )
    click.echo(datetime.datetime.now().strftime('End time: %Y-%m-%d %H:%M'))


@cli.command(cls=click.RichCommand, context_settings={'show_default': True})
@click.argument(
    'directory', type=click.Path(exists=True, file_okay=False, dir_okay=True)
)
@click.option(
    '--ref-labels',
    type=str,
    default='noa',
    help='Reference labels.',
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
    help='Output file for the report.',
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


@cli.command(cls=click.RichCommand)
@click.argument(
    'files_or_dir',
    nargs=-1,
    type=click.Path(exists=True),
    callback=_files_or_dir,
)
@click.option(
    '--codons/--no-codons',
    '-c/-C',
    default=False,
    help='Group the nucleotides by codon',
    show_default=True,
)
@click.option(
    '--aa/--no-aa',
    '-a/-A',
    default=False,
    help='Display the translated peptide',
    show_default=True,
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
