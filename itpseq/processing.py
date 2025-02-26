import re
import json
import pandas as pd
import numpy as np
from pathlib import Path

from functools import lru_cache
from typing import overload, Optional, Union, Any
from pandas import DataFrame

from .utils import fcache

__all__ = [
    'read_aafile_as_series',
    'read_log_json',
    'code2pos',
    'ranges',
    'get_ribosome_site',
    'compute_counts',
    'itp_len_add_positions',
    'DE',
    'flag_low_density',
]

# @overload
@lru_cache
@fcache
# @log
def read_aafile_as_series(
    filename,
    min_peptide=None,
    max_peptide=None,
    how='aax',
    limit=None,
    sample=None,
):
    """Reads the aminoacid inverse-toeprint file as a pandas Series, filters entries based on peptide length and stop codons.

    Parameters
    ----------
    filename : Path or str
        Path to the amino-acid file
    min_peptide : int, optional
        Minimum peptide length to keep, by default None
    max_peptide : int, optional
        Maximum peptide length to keep, by default None
    how : str, optional
        Mode to filter the stops: "aax" will remove peptides with a stop before the A-site, by default 'aax'

    Returns
    -------
    Series
        _description_
    """
    with open(filename) as f:
        if min_peptide or max_peptide or how == 'aax':
            S = '*' if how == 'aax' else ''
            min_ = min_peptide - 1 if min_peptide else ''
            max_ = max_peptide - 1 if max_peptide else ''
            pat = re.compile(rf'^\s*[^\s{S}]{{{min_},{max_}}}.$')
            out = pd.Series(
                [
                    line
                    for line in (l.rstrip('\n') for l in f.readlines())
                    if not line.startswith('#') and pat.match(line)
                ]
            )
        else:
            out = pd.Series(
                [
                    line.rstrip('\n')
                    for line in f.readlines()
                    if not line.startswith('#')
                ]
            )

        if limit:
            out = out.head(limit)
        if sample:
            out = out.sample(sample)

        return out


def read_log_json(path: Union[str, Path]):
    """Reads the metadata from the fastq parsing for this replicate.

    Parameters
    ----------
    path : Union[str,Path]
        Path to the metadata json file.

    Returns
    -------
    data: dict
        Miscellaneous information on the sequences.
    """

    data = json.load(open(path))
    return data

    # df_data = {}
    # for k in list(data):
    #     if isinstance(data[k], dict):
    #         df_data[k] = data.pop(k)

    # return (data, pd.DataFrame(df_data)
    #                 .rename(int)
    #                 .sort_index()
    #                 .convert_dtypes()
    #        )


def code2pos(code):
    """
    Converts a ribosome site or code to the python on the peptide relative to the end (A-site = -1).

    Parameters
    ----------
    code : str or int
        Ribosome site (E, P, A) or position relative to the P-site.

    Returns
    -------
    int
        Matching position relative to the end of the peptide (A-site = -1, P-site = -2)
    """
    d = {'E': -1, 'P': 0, 'A': +1}
    code = d.get(code, code)
    try:
        return (
            int(code) - 2
        )  # A-site position is +1, python relative to end = -1
    except ValueError:
        raise ValueError(f'{code} is an invalid code')


def ranges(s=''):
    """transforms 'x-y' into slice(x, y, None)"""
    return [
        slice(y[0], y[1] + 1 if y[1] != -1 else None)
        if len(y := list(map(code2pos, x.split(':')))) > 1
        else y[0]
        for x in s.split(',')
    ]


def get_ribosome_site(pos, short=False):
    """
    Maps a position relative to the ribosome to its corresponding ribosomal site.

    This function translates a given position (e.g., `-1`, `0`, `1`, or `E`, `P`, `A`)
    into its corresponding ribosomal site (E-site, P-site, A-site) or a numerical representation
    if the position does not correspond to a standard site.

    Parameters
    ----------
    pos : int or str
        The position relative to the ribosome. Can be:
        - `-1`, `0`, `1` for standard ribosomal sites.
        - `E`, `P`, `A` for explicit site designations.
        - Other numeric values, which will be returned in `±<absolute value>` format.
    short : bool, optional
        If True, returns shorter labels (`E`, `P`, `A`). Defaults to False,
        which returns full names (e.g., `E-site`, `P-site`, `A-site`).

    Returns
    -------
    str
        The ribosomal site as a string. If the input does not match a standard site,
        the function returns the input position formatted as `±<absolute value>`.
        In case of invalid input, the same value is returned.

    Notes
    -----
    - Valid inputs include both numeric (`-1`, `0`, `1`) and alphabetical ribosomal site labels (`E`, `P`, `A`).
    - For numeric inputs outside the standard range, the output is formatted as `+<value>` or `−<value>`.
    """

    pos = str(pos)

    # fmt: off
    positions = {'-1': 'E', '0': 'P', '1': 'A',
                 'E': 'E', 'P': 'P', 'A': 'A',
                 }
    # fmt: on
    if not short:
        positions = {k: f'{v}-site' for k, v in positions.items()}
    try:
        if pos in positions:
            return positions[str(pos)]
        else:
            return f'{"+" if (pos := int(pos)) > 0 else "−"}{abs(pos)}'
    except TypeError:
        print(f'Error in get_ribosome_site({pos})')
        return pos


@fcache
def compute_counts(seqs_series, pos):
    if not pos:
        return (
            seqs_series.str.split('(?<=.)(?=.)', expand=True)
            .apply(lambda x: x.value_counts())
            .pipe(lambda x: x.rename(columns=lambda c: c - x.shape[1] + 2))
            # .pipe(lambda x: x.set_axis(np.arange(x.shape[1])-x.shape[1]+2, axis=1))
        )

    from functools import reduce

    slices = ranges(pos)
    return reduce(
        lambda x, y: x + y, map(lambda x: seqs_series.str[x], slices)
    ).value_counts()


def itp_len_add_positions(ax, min_codon=0, max_codon=10):
    greys = ['#808080', '#CACACA']
    start = 18
    ax.axvline(start + 3 * min_codon, color='k', alpha=0.1, zorder=-1)
    for i in range(min_codon, max_codon):
        ax.axvspan(
            (start + 3 * i),
            (start + 3 * i + 3),
            color=greys[i % 2] if i else '#D4C13A',
            ls='',
            alpha=0.1,
            zorder=-1,
        )
        ax.axvline((start + 3 * i + 3), color='k', alpha=0.1, zorder=-1)
        ax.text(
            start + 3 * i + 1.5,
            ax.get_ylim()[1] * 0.99,
            i + 1,
            ha='center',
            va='top',
            zorder=-1,
            color=greys[0],
        )
    return ax


@fcache
def DE(
    df: DataFrame,
    sample_df: DataFrame,
    cond: str = None,
    ref: str = None,
    join: bool = False,
    quiet: bool = True,
    multi: bool = True,
    n_cpus: int = None,
    raw: bool = False,
) -> DataFrame | None | Any:
    """
    Performs differential expression (DE) analysis using DESeq2 via the `pydeseq2` library.

    This function calculates enrichment statistics between a condition and a reference sample,
    leveraging the replicates of both samples.

    Parameters

    ----------
    df : pandas.DataFrame
        A DataFrame containing counts data, where rows represent the positions/motifs
        and columns correspond to samples.
    sample_df : pandas.DataFrame
        A DataFrame providing metadata for the samples.
        Must include a column with sample names matching `df` and a `sample` column (used as the design factor).
    cond : str
        The name of the condition to analyze (used in contrast).
    ref : str
        The reference condition to compare against (used in contrast).
    join : bool, optional
        If True, joins the DE results back to the original `df`. Defaults to False.
    quiet : bool, optional
        If True, suppresses the console output of the `pydeseq2` library. Defaults to True.
    multi : bool, optional
        Whether to compute DE with a specific contrast (`cond` vs. `ref`). Defaults to True.
    n_cpus : int, optional
        The number of CPUs to utilize for parallel processing. Defaults to the total number of available CPUs.

    Returns
    -------
    pandas.DataFrame or None or Any
        - If successful and `join` is False, returns a DataFrame containing DE results, sorted by the `log2FoldChange` column.
        - If `join` is True, returns the input `df` with the DE results joined as additional columns.
        - If the analysis fails, returns None.

    Notes
    -----
    - The DE results include additional columns including:
      - `log2FoldChange`: Log2 fold change enrichment of the position/motif relative to the reference.
      - `pvalue`, `padj`: P-values and adjusted P-values.
      - `log10pvalue`, `log10padj`: -Log-transformed P-values for visualization.

    Exceptions
    ----------
    - Prints a message and returns `None` if the DESeq2 analysis fails due to invalid input or other issues.
    """

    import contextlib, io

    f = io.StringIO()

    if not n_cpus:
        import multiprocessing

        n_cpus = multiprocessing.cpu_count()

    cms = []
    if quiet:  # if quiet, catch the stdout/stderr of pydeseq2
        f = io.StringIO()
        cms = [contextlib.redirect_stdout(f), contextlib.redirect_stderr(f)]

    with contextlib.ExitStack() as es:
        for cm in cms:
            es.enter_context(cm)

        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        dds = DeseqDataSet(
            counts=df.fillna(0, downcast='infer')[sample_df['samplename']].T,
            metadata=sample_df[['sample']],
            design_factors='sample',
            refit_cooks=True,
            n_cpus=n_cpus,
        )

        try:
            dds.deseq2()

            if multi:
                stat_res = DeseqStats(dds, contrast=['sample', cond, ref])
            else:
                stat_res = DeseqStats(dds)
            stat_res.summary()

            if raw:
                return dds, stat_res

            res = stat_res.results_df.sort_values(
                by='log2FoldChange', ascending=False
            )

            res['log10pvalue'] = -np.log10(res['pvalue'])
            res['log10padj'] = -np.log10(res['padj'])

            if join:
                return df.join(res, how='right')
            return res
        except ValueError:
            print(f'DESeq2 failed for {cond}')
            return None


def flag_low_density(x, y, bins=100, thresh=None):
    """Counts the points on a 2D grid and flags those in bins with low counts"""
    hist, bins_x, bins_y = np.histogram2d(x, y, bins=bins)
    bx = np.digitize(x, bins=bins_x, right=True) - 1
    by = np.digitize(y, bins=bins_y, right=True) - 1
    if thresh:
        return hist[bx, by] <= thresh
    return hist[bx, by]
