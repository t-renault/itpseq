"""Helper functions to produce graphs from the itpseq data"""

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.collections import LineCollection

from .utils import aa_colors

__all__ = ['itoeprint_plot', 'motif_logo']


def itoeprint_plot(
    dataset=None,
    plot='bands',
    norm='sum',
    norm_range=(21, 51),
    exposure=1,
    limit=(0, 100),
    show_range=True,
    ax=None,
):
    """
    Plots a virtual inverse-toeprint gel.

        Parameters
        ----------
        dataset : str, optional
            DataSet or Sample to use as source of data (needs a ``toeprint_df`` attribute)
        plot : str, optional
            Type of representation. "bands" will display lines with a width proportional to the number of reads.
            "shades" will display a fixed-height line with a varying shade.
        norm : str, optional
            Type of normalization for each lane (possible values are None, 'mean', 'median', 'max', 'std')
        norm_range : tuple, optional
            Range of lengths considered to perform the normalization.
            Defaults to 21-51, which corresponds to the first 10 codons after the start.
        exposure : int, optional
            Modulates the global intensity of the bands in the "bands" type of plot.
        limit: tuple, optional
            Limits the visible range of lengths.
        show_range : bool, optional
            Shows the regions excluded from the normalization in light red.
        ax : matplotlib.axes.Axes, optional
             ax to use for plotting, otherwise create a new figure.
    """
    if not ax:
        ax = plt.subplot()

    df = dataset.toeprint_df

    start, end = norm_range
    low_limit, high_limit = limit

    df = df.div(df.loc[end:start].agg(norm) if norm else df.max().max() / 20)

    xticks_kwargs = dict(
        rotation='vertical' if df.shape[1] > 10 else 'horizontal', ha='center'
    )

    if plot == 'bands':
        for i, c in enumerate(df):
            s = df[c].dropna()
            points = [[(i - 0.4, y), (i + 0.4, y)] for y in s.index]
            # ax.add_collection(LineCollection(points, color='k', lw=s*10, alpha=0.3))
            ax.add_collection(
                LineCollection(
                    points,
                    color='k',
                    lw=s * 0.2 * exposure,
                    path_effects=[path_effects.Stroke(capstyle='round')],
                    # path_effects=[path_effects.Stroke(capstyle='butt')],
                )
            )

        ax.set_xlim(-0.5, df.shape[1] - 0.5)
        ax.set_ylim(low_limit, high_limit)
        ax.set_xticks(np.arange(df.shape[1]), df.columns, **xticks_kwargs)

        ax.spines[['left', 'right', 'top', 'bottom']].set_visible(False)

        if show_range:
            ax.axhline(start, color='red', alpha=0.5)
            ax.axhline(end, color='red', alpha=0.5)
            ax.axhspan(low_limit, start, color='red', alpha=0.1, zorder=2)
            ax.axhspan(end, high_limit, color='red', alpha=0.1, zorder=2)

    elif plot == 'shades':
        sns.heatmap(
            df.pipe(
                lambda x: x.reindex(
                    range(high_limit or df.index.max(), low_limit, -1)
                )
            ),
            cmap='Greys',
            ax=ax,
            cbar=False,
        )
        i = 0   # avoid static checker complaint
        for i in range(0, df.shape[1] + 1):
            ax.axvline(i, c='w', lw=70 / (df.shape[1] + 1), zorder=1)

        ax.set_xlim(0, i)
        ax.set_xticks(
            np.arange(df.shape[1]) + 0.5, df.columns, **xticks_kwargs
        )

        if show_range:
            ax.axhline(high_limit - start, color='red', alpha=0.5)
            ax.axhline(high_limit - end, color='red', alpha=0.5)
            ax.axhspan(
                high_limit,
                high_limit - start,
                color='red',
                alpha=0.1,
                zorder=1,
            )
            ax.axhspan(
                high_limit - end, low_limit, color='red', alpha=0.1, zorder=1
            )

    ax.set_ylabel('Inverse-toeprint length')
    ax.tick_params(top=False, labeltop=True, bottom=False)

    return ax


def motif_logo(
    df,
    ref_spl=None,
    *,
    ax=None,
    logo_type='extra_counts_bits',
    logo_kwargs={'color_scheme': aa_colors},
    return_matrix=False,
):
    """Plots a logo based on a DatFrame of the sample/reference counts.

    Parameters
    ----------
    df: pandas.DataFrame
        The input DataFrame must have the motifs as index and two columns "sample" and "ref" which contain the
        average normalized counts of the sample and reference.
    ref_spl: list or None,
        Instead of passing a DataFrame with "sample" and "ref" columns, it is also possible to define the column names
         to use by passing a list of the [ref, sample] column names.
    ax: matplotlib Axes, optional
        If passed, the figure will be drawn on the given axes.
    logo_type: str, optional
        Type of logo to compute:
        * "raw_freq": unweighted frequencies of the amino-acids for all present motifs
        * "extra_counts": Computes the sum of extra counts (sample - reference) for each residue per position
        * "sum_log2FC": sum of the log2FoldChange for each residue per position
        * "<logo_type>_bits": If any of the above has a "_bits" suffix, an extra conversion to bits is performed.

    return_matrix: bool, optional
        If True, the logo matrix is returned as together with the logo as (logo, matrix).
    """
    if df.empty:
        print('Input DataFrame is empty. Cannot compute a logo.')
        return
    assert logo_type in {
        'raw_freq',
        'extra_counts',
        'sum_log2FC',
        'raw_freq_bits',
        'extra_counts_bits',
        'sum_log2FC_bits',
    }

    import logomaker

    if isinstance(ref_spl, list) and len(ref_spl) == 2:
        # df = df[ref_spl].set_axis(['ref', 'sample'], axis='columns')
        df = df.rename(columns=dict(zip(ref_spl, ['ref', 'sample'])))

    df = df.rename_axis('motif')

    if not ax:
        ax = plt.subplot()

    if logo_type in {'raw_freq', 'raw_freq_bits'}:
        # counts the frequency of each amino-acid at each position
        # no weighting is done based on the enrichment data
        tmp = df.index.to_series().str.extractall('(?P<aa>.)').reset_index()
        matrix = pd.crosstab(tmp['match'], tmp['aa'], normalize='index')
        if logo_type == 'raw_freq_bits':
            matrix = logomaker.transform_matrix(
                matrix, from_type='probability', to_type='information'
            )
    elif logo_type in {'extra_counts', 'extra_counts_bits'}:
        # computes the sum of extra counts (sample - ref)
        w = df['sample'].sub(df['ref']).rename('weight')
        cnts = (
            df.index.to_series()
            .str.extractall('(?P<aa>.)')
            .join(w)
            .reset_index()
            .pivot_table(
                index='match', columns='aa', values='weight', aggfunc='sum'
            )
        )
        mask = cnts.isna()
        matrix = cnts.fillna(0)
        if logo_type == 'extra_counts':
            matrix = logomaker.transform_matrix(
                matrix, from_type='counts', to_type='probability'
            )
        elif logo_type == 'extra_counts_bits':
            matrix = logomaker.transform_matrix(
                matrix, from_type='counts', to_type='information'
            )
        # logomaker transforms zeros into small numbers
        matrix = matrix.mask(mask, 0)
    # elif logo_type in {'sum_log2FC', 'sum_log2FC_bits'}:
    else:
        matrix = (
            df.index.to_series()
            .str.extractall('(?P<aa>.)')
            .join(df['log2FoldChange'])
            .reset_index()
            .pivot_table(
                index='match',
                columns='aa',
                values='log2FoldChange',
                aggfunc='sum',
            )
        ).fillna(0)
        matrix = matrix.div(matrix.sum(axis=1), axis=0)
        if logo_type == 'sum_log2FC_bits':
            matrix = logomaker.transform_matrix(
                matrix, from_type='probability', to_type='information'
            )

    logo = logomaker.Logo(matrix, **logo_kwargs, ax=ax)
    if return_matrix:
        return logo, matrix
    return logo
