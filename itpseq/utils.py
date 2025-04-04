"""General helper functions for the itpseq library"""

import io
import logging
import tempfile
from collections import defaultdict
from functools import wraps
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from pandas.api.types import is_numeric_dtype
from seaborn.axisgrid import FacetGrid

__all__ = [
    'log',
    'fcache',
    'aa_colors',
    'aa_order',
    'plot_to_html',
    'table_to_html',
    'dict_to_tuple',
    'dedup_names',
    'itp_schematic',
]

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    filename=Path(tempfile.gettempdir()) / 'itpseq.log',
    filemode='a',
)
logger = logging.getLogger('itpseq')
logger.addHandler(logging.StreamHandler())


def dict_to_tuple(d, *, ignore=None, keep=None):
    """Converts a dictionary to a tuple of sorted keys and values"""
    if ignore or keep:
        if isinstance(ignore, str):
            ignore = {ignore}
        if keep is not None and isinstance(keep, str):
            keep = {keep}
        return tuple(
            sorted(
                (k, v)
                for k, v in d.items()
                if (ignore is None or k not in ignore)
                and (keep is None or k in keep)
            )
        )
    return tuple(sorted(d.items()))


def dedup_names(names, sep='.'):
    """
    Rename labels if duplicates exist.

    Inspired by pandas.io.common.dedup_names

    Parameters
    ----------
    names :
        Iterable of names to deduplicate.

    sep : str, optional
        String to use a separator with the suffix.

    Examples
    --------
    >>> dedup_names(['x', 'y', 'x', 'x'])
    ['x', 'y', 'x.1', 'x.2']
    """
    names = list(names)
    counts = defaultdict(int)

    for i, label in enumerate(names):
        cur_count = counts[label]
        while cur_count > 0:
            counts[label] = cur_count + 1
            label = f'{label}{sep}{cur_count}'
            cur_count = counts[label]
        names[i] = label
        counts[label] = cur_count + 1

    return names


def log(func):
    """Decorator to log function calls"""

    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.info(f'{func.__name__}({args}, {kwargs})')
        return func(*args, **kwargs)

    return wrapper


def fcache(func):
    """Decorator to cache function calls returning DataFrames to csv files"""

    @wraps(func)
    def wrapper(
        *args,
        _cache_dir=None,
        _cache_prefix=None,
        _nocache=False,
        _force=False,
        **kwargs,
    ):
        if _nocache or _cache_dir is None:
            return func(*args, **kwargs)
        Path(_cache_dir).mkdir(parents=True, exist_ok=True)
        # args_str = '' '_' + '_'.join(map(str, args)) if args else ''
        args_str = ''   # FIXME decide how to handle filenames as parameters
        kwargs_str = '_' + '_'.join(
            f'{k}={v}' for k, v in sorted(kwargs.items())
        )
        _cache_prefix = f'_{_cache_prefix}' if _cache_prefix else ''
        cache_file = (
            Path(_cache_dir)
            / f'{func.__name__}{_cache_prefix}{args_str}{kwargs_str}.csv'
        )
        if (not _nocache) and cache_file.exists():
            print(f'Loading cached {cache_file}')
            result = pd.read_csv(
                cache_file, index_col=0, na_values=[''], keep_default_na=False
            )
            if result.shape[1] == 1 and result.columns[0] in {'0', 'count'}:
                result = result.squeeze(axis=1).rename(
                    _cache_prefix
                )  # we want a Series
            return result
        result = func(*args, **kwargs)
        if result is not None:
            result.to_csv(cache_file)
        return result

    return wrapper


# amino acids color codes
HYDRO = '#74C170'      # hydrophobic
SPECIAL = '#E4DF51'    # special
NEGATIVE = '#7094C1'   # negatively charged
POSITIVE = '#DB4755'   # positively charged
POLAR = '#AF70C1'      # polar

aa_colors = {
    'A': HYDRO,
    'C': SPECIAL,
    'D': NEGATIVE,
    'E': NEGATIVE,
    'F': HYDRO,
    'G': SPECIAL,
    'H': POSITIVE,
    'I': HYDRO,
    'K': POSITIVE,
    'L': HYDRO,
    'M': HYDRO,
    'N': POLAR,
    'P': SPECIAL,
    'Q': POLAR,
    'R': POSITIVE,
    'S': POLAR,
    'T': POLAR,
    'V': HYDRO,
    'W': HYDRO,
    'Y': HYDRO,
    '*': 'grey',
    'm': HYDRO,
}

del HYDRO, SPECIAL, NEGATIVE, POSITIVE, POLAR

aa_order = list('HRKDESTNQCGPAVILMFYW*m')


def plot_to_html(
    fig, format='svg', figsize=None, dpi=150, destroy=True, bgcolor='none'
):
    """Converts a matplotlib figure or axes to an SVG or PNG image that can be
    embedded in HTML"""

    if format not in ('svg', 'png'):
        raise ValueError(f'format must be one of "svg" or "png", got {format}')

    if isinstance(fig, Figure):
        pass
    elif isinstance(fig, Axes):
        fig = fig.get_figure()
    elif isinstance(fig, FacetGrid):
        fig = fig.fig
    elif fig is None:
        return None
    else:
        return f'unknown object type: {type(fig)}'

    facecolor = None
    if bgcolor:
        facecolor = fig.patch.get_facecolor()
        fig.patch.set_facecolor(bgcolor)

    if figsize:
        fig.set_size_inches(*figsize)

    if format == 'svg':
        out = io.StringIO()
        fig.savefig(out, format='svg', bbox_inches='tight', pad_inches=0.1)
        outstr = out.getvalue()
    elif format == 'png':
        import base64

        out = io.BytesIO()
        fig.savefig(
            out, format='png', dpi=dpi, bbox_inches='tight', pad_inches=0.1
        )
        b64 = base64.b64encode(out.getvalue()).decode('utf-8')
        outstr = f"<img src='data:image/png;base64,{b64}'>"
    else:
        return None

    if destroy:
        plt.close(fig)
    else:
        if facecolor:
            fig.patch.set_facecolor(
                facecolor
            )   # restore original background color
    return outstr


def table_to_html(df):
    """Converts a Pandas DataFrame to a HTML table"""
    return (
        df.style.set_properties(**{'text-align': 'right'})
        .set_table_styles(
            [{'selector': 'th', 'props': [('text-align', 'right')]}]
        )
        .format(
            {
                c: lambda x: f'{x:,.0f}'.replace(',', ' ')
                for c in df
                if is_numeric_dtype(df[c]) and df[c].max() > 1_000
            }
        )
        .to_html()
    )


def itp_schematic(
    *, codons=9, coding_sequence=None, min_overhang=12, width_factor=1
):
    """
    Creates a schematic of an iTP-Seq read with annotations and positions.

    Parameters
    ----------
    codons : int
        Number of codons to use for the translated sequence.
        This is ignored if an explicit coding sequence is provided.
    coding_sequence : str
        Explicit coding sequence to use.
    min_overhang : int
        Minimum length of the overhang (up to two extra nucleotides can be present).
    width_factor : float
        Factor to change the schematic width.

    Returns
    -------
    matplotlib.Axes

    Examples
    --------

    Produce a default example with a dummy sequence:
     >>> itp_schematic()

    .. image:: /_static/itp_schematic.png

    Produce a schematic with only 2 codons:
     >>> itp_schematic(codons=2)

    .. image:: /_static/itp_schematic_2codons.png

    Produce a schematic with a specific sequence:
     >>> itp_schematic(sequence='ATGGCACCTCTAGAG')

     .. image:: /_static/itp_schematic_custom.png
    """
    import numpy as np
    from dna_features_viewer import GraphicFeature, GraphicRecord
    from matplotlib.ticker import FuncFormatter

    extra = 2  # number of possible extra nucleotides
    width_factor *= 0.3

    if coding_sequence:
        assert (
            len(coding_sequence) % 3 == 0
        ), f'coding_sequence should have a multiple of 3 nucleotides, got {len(coding_sequence)}'

        sequence = coding_sequence + 'n' * min_overhang + '-' * extra
        codons = len(coding_sequence) // 3
        translation = None
    else:
        sequence = 'N' * 3 * codons + 'n' * min_overhang + '-' * extra
        translation = 'm' * (codons > 0) + 'X' * (codons - 1)

    from_P = 6 + min_overhang + extra  # number of nt after from P-site (0)
    offset = len(sequence) - from_P

    labels = ['E', 'P', 'A']

    record = GraphicRecord(
        sequence=sequence,
        features=[
            GraphicFeature(
                start=(codons - len(labels) + i) * 3,
                end=(codons - len(labels) + 1 + i) * 3,
                fontdict={'weight': 'bold'},
                label=labels[i],
                strand=0,
                color='#ffcccc',
            )
            for i in range(len(labels))
        ]
        + [
            GraphicFeature(
                start=len(sequence) - min_overhang - extra,
                end=len(sequence),
                strand=0,
                color='w',
                label='protected overhang',
            )
        ],
        labels_spacing=0,
    )

    # +9 is an offset to account for the extra space on the side of the figure
    ax, _ = record.plot(
        figure_width=width_factor * (len(sequence) / 3 + 9), figure_height=1.5
    )
    record.plot_sequence(ax)
    record.plot_translation(
        ax,
        (0, codons * 3),
        y_offset=-1,
        translation=translation,
        long_form_translation=False,
    )

    xticks = (
        np.array(
            [
                *range(-(codons - 2) * 3, -6, 3),
                *range(-3 if codons < 4 else -6, 6),
                *range(3, from_P, 3),
            ]
        )
        + offset
    )
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: x - offset))
    ax.set_xticks(xticks)
    ax.set_ylim(-1, 1)

    ax2 = ax.secondary_xaxis('top')
    ax2.set_xticks(range(-8, codons * 3, 3))
    ax2.xaxis.set_major_formatter(
        FuncFormatter(lambda x, pos: (x + 2) // 3 - codons + 1)
    )
    ax2.tick_params(length=0, colors='grey')
    ax2.spines[['top']].set_visible(False)

    if codons < 3:
        ax.set_xlim(left=ax.get_xlim()[0] - 0.8)

    return ax
