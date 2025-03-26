"""Entry point for the command line tool"""

import sys
import webbrowser
from pathlib import Path

import click

from . import __version__, parsing
from .core import DataSet


@click.group(
    context_settings={'help_option_names': ['-h', '--help']},
    epilog='To get help on a specific command, run: itpseq COMMAND --help',
)
@click.version_option(__version__, '--version', '-v', prog_name='itpseq')
def main():
    """itpseq - A command-line tool for the analysis of iTP-Seq data."""
    print(
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


@main.command()
def help():
    """Open the HTML documentation in the default web browser."""

    module_dir = Path(__file__).resolve().parent.parent
    docs_index = module_dir.joinpath('docs', '_build', 'html', 'index.html')

    if not docs_index.exists():
        print(
            'Error: Documentation not found. Please build the Sphinx documentation first.'
        )
        return

    # Open the documentation in the default web browser
    webbrowser.open(f'file://{docs_index}')


@main.command()
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


@main.command(
    context_settings={'ignore_unknown_options': True}, add_help_option=False
)
@click.argument('args', nargs=-1)
def parse(args):
    """
    Parse and filter the assembled iTP-Seq FASTQ files to produce files needed for the analysis.
    """
    # Pass arguments directly to parsing.main()
    sys.argv = ['parse'] + list(args)  # Simulate sys.argv for argparse
    parsing.main()


if __name__ == '__main__':
    main()
