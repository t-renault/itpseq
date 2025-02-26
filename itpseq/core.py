from pathlib import Path
from typing import Optional, Union
from collections import defaultdict

from functools import total_ordering, lru_cache, cached_property

import re
import datetime
from copy import copy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from .utils import *
from .processing import *
from .plotting import *

IPYTHON = False
try:
    from IPython import get_ipython

    ip = get_ipython()
    if ip is not None and getattr(ip, 'kernel', None) is not None:
        IPYTHON = True
except ImportError:
    pass

__all__ = ['DataSet', 'Sample', 'Replicate']


@total_ordering
class Replicate:
    """
    Replicate instances represent a specific biological or experimental replicate of a Sample.

    The purpose of the class is to handle, process, and analyze data corresponding to a replicate.
    Replicate objects provide methods to load associated data, compute statistical measures,
    and generate graphical representations such as sequence logos.

    Attributes
    ----------
    filename : Optional[Path]
        Path to the file associated with the replicate. This file is expected to
        contain raw data relevant to the replicate.
    sample : Optional[Sample]
        The sample object this replicate belongs to.
    replicate : Optional[str]
        Identifier or label for the replicate (e.g., "1").
    labels : Optional[dict]
        Dictionary of labels or metadata associated with the replicate.
    name : str
        Name of the sample, derived from `sample.name` if provided.
    dataset : Any
        The DataSet the Sample belongs to, derived from `sample.dataset` if provided.
    kwargs : dict
        Additional keyword arguments and metadata stored as "meta" during initialization.
    """

    def __init__(
        self,
        *,
        replicate: Optional[str] = None,
        filename: Optional[Path] = None,
        sample: Optional['Sample'] = None,
        labels: Optional[dict] = None,
        **kwargs,
    ):
        # print(f'Replicate {kwargs=}')
        self.filename = filename
        # self.file_pattern = file_pattern or '{lib_type}_{sample}{replicate}_aa.processed.txt'
        self.sample = sample
        self.replicate = replicate
        self.labels = labels
        self.sample_name = self.sample.name if sample else ''
        self.name = (
            f'{self.sample_name}.{self.replicate}'
            if self.sample
            else self.replicate
        )
        self.dataset = self.sample.dataset if sample else None
        self._cache_dir = self.dataset.cache_path if self.dataset else None
        self.meta = kwargs
        #self._raw_data = None

    def __repr__(self):
        return f'Replicate({self.name})'

    def __str__(self):
        return self.name

    def __gt__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return (self.sample, self.replicate) > (other.sample, other.replicate)

    def load_data(self, **kwargs):
        return read_aafile_as_series(
            self.filename,
            _cache_dir=self._cache_dir,
            _cache_prefix=self.name,
            **kwargs,
        )

    @lru_cache
    def get_counts(self, pos=None, **kwargs):
        prefix = (
            '_' + '_'.join(f'{k}={v}' for k, v in sorted(kwargs.items()))
            if kwargs
            else ''
        )
        return compute_counts(
            self.load_data(**kwargs),
            pos=pos,
            _cache_dir=self._cache_dir,
            _cache_prefix=self.name + prefix,
        )

    def logo(
        self,
        logo_kwargs=None,
        ax=None,
        fMet=False,
        type='information',
        **kwargs,
    ):
        """
        Generates a sequence logo based on the aligned inverse-toeprints, using the `logomaker` library.

        Parameters
        ----------
        logo_kwargs : dict, optional
            Additional keyword arguments passed to `logomaker.Logo` for customizing
            the sequence logo. Defaults to `{'color_scheme': 'NajafabadiEtAl2017'}`.
        ax : matplotlib.axes.Axes, optional
            Pre-existing matplotlib Axes to draw the logo on. A new Axes is created if not provided.
        fMet : bool, optional
            If False, removes `m` (formyl-methionine / start codon) from the alignment when building the logo. Defaults to False.
        type : str, optional
            The transformation type applied to the counts matrix. Possible values include:
            - `'information'` for information content.
            - `'probability'` for probabilities.
            Defaults to `'information'`.
        **kwargs : dict
            Additional keyword arguments passed to filter the input data (e.g., `pos`, `min_peptide`, `max_peptide`...).

        Returns
        -------
        logomaker.Logo
            A `logomaker.Logo` object representing the sequence logo.

        Notes
        -----
        - Sequence alignment data is first converted to a counts matrix via the `logomaker.alignment_to_matrix` method.
        - The ribosomal site corresponding to each position is annotated on the x-axis.
        - Transformation of the counts matrix (e.g., `counts` to `information`) is performed using `logomaker.transform_matrix`.

        Examples
        --------
        # Simple logo plot with default settings
        logo = obj.logo()

        # Logo plot with min_peptide filtering
        logo = obj.logo(min_peptide=3)

        # Logo plot with custom transformation type and filtering
        logo = obj.logo(type='probability', min_peptides=2, fMet=True)
        """

        if logo_kwargs is None:
            logo_kwargs = {'color_scheme': 'NajafabadiEtAl2017'}
        import logomaker

        df = (
            logomaker.alignment_to_matrix(self.load_data(**kwargs))
            # .pipe(lambda x: x.rename(lambda c: get_ribosome_site(c-x.shape[0]+2)))
            .drop(columns=' ', errors='ignore')
        )
        if not fMet:
            df = df.drop(columns='m', errors='ignore')

        df = logomaker.transform_matrix(df, from_type='counts', to_type=type)
        # print(df)

        logo = logomaker.Logo(df, ax=ax, **logo_kwargs)
        logo.ax.set_xticks(
            range(df.shape[0]),
            (df.index - df.shape[0] + 2).map(get_ribosome_site),
        )
        return logo


@total_ordering
class Sample:
    """
    Represents a sample in a dataset, its replicates, reference, and associated metadata.

    The `Sample` class is used to encapsulate information and behavior related to samples in a dataset.
    It manages details like labels, references, replicates, and metadata, and provides methods for analyzing
    replicates, performing differential enrichment analysis, and creating visualizations.
    """

    def __init__(
        self,
        *,
        labels: dict,
        reference=None,
        dataset=None,
        data=None,
        keys=('sample',),
        **kwargs,
    ):
        """ """
        # print(f'Sample {kwargs=}')
        self.labels = labels
        self.name = '.'.join(labels[k] for k in keys)
        self.dataset = dataset
        self._cache_dir = self.dataset.cache_path if self.dataset else None
        self.reference = reference  # If None, this sample is the reference
        self.data = data
        self.meta = kwargs
        # self.replicates = []
        #
        # for p, rep in self.paths:
        #    self.replicates.append(Replicate(filename=p, replicate=rep, sample=self))

        # self.replicates = sorted([Replicate(filename=p, replicate=rep, sample=self, labels=labels)
        #                          for p, rep, labels in self.paths])
        self.replicates = sorted(
            [Replicate(**data, sample=self) for data in self.data]
        )

    def __getitem__(self, key):
        return self.replicates[key]

    def __repr__(self):
        r = (
            f':[{", ".join([r.replicate for r in self.replicates])}]'
            if self.replicates
            else ''
        )
        ref = f', ref: {self.reference}' if self.reference else ''
        return f'Sample({self.name}{r}{ref})'

    def __str__(self):
        return self.name

    def __truediv__(self, other):
        return MetaSample(sample=self, reference=other)

    def __gt__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        if other == self.reference:
            return True
        elif other.reference == self:
            return False
        return (self.dataset, self.name) > (other.dataset, other.name)

    @property
    def name_vs_ref(self):
        return (
            f'{self.name} vs {self.reference.name}'
            if self.reference
            else self.name
        )

    @property
    def name_ref(self):
        return (
            f'{self.name}-{self.reference.name}'
            if self.reference
            else self.name
        )

    def load_replicates(self): ## FIXME is this useful to keep?
        for replicate in self.replicates:
            replicate.load_data()

    def infos(self, html=False):
        order = [
            'total_sequences',
            'noadaptor',
            'contaminant',
            'lowqual',
            'tooshort',
            'toolong',
            'extra0',
            'extra1',
            'extra2',
        ]
        key = {c: i for i, c in enumerate(order)}
        out = (
            pd.DataFrame.from_dict(
                {r.name: r.meta['fastq'] for r in self.replicates},
                orient='index',
            )
            .sort_index(axis=1, key=lambda x: x.map(key))
            .select_dtypes('number')
        )
        if html:
            return table_to_html(out)
        return out

    @lru_cache
    def get_counts(self, pos=None, **kwargs):
        # print(locals())
        return pd.concat(
            {str(r): r.get_counts(pos, **kwargs) for r in self.replicates},
            axis=1,
        )

    # @lru_cache
    def get_counts_ratio(
        self, pos=None, factor=1_000_000, exclude_empty=True, **kwargs
    ):
        counts = self.get_counts(pos=pos, **kwargs)  # compute position counts
        if exclude_empty:
            counts = counts[~counts.index.str.fullmatch(' *')]
        norm = counts.div(counts.sum()).mul(factor)  # normalize
        avg = norm.mean(axis=1)
        out = counts.assign(**{self.name: avg})
        if self.reference:
            ref = (
                self.reference.get_counts_ratio(
                    pos=pos, factor=factor, **kwargs
                ).filter(like=self.reference.name)
                # .drop(columns='ratio', errors='ignore')
            )
            cols = [
                *ref.columns[:-1],
                *out.columns[:-1],
                *ref.columns[-1:],
                *out.columns[-1:],
            ]

            # print(ref, out)

            return ref.join(out)[cols].assign(
                ratio=out[self.name] / ref[self.reference.name]
            )
        else:
            return out
            # normalize

    @lru_cache
    def DE(
        self,
        pos=None,
        join=False,
        quiet=True,
        filter_size=True,
        multi=True,
        n_cpus=None,
        raw=False,
        _nocache=False,
        **kwargs,
    ):
        print(f'DE: {self.name}, {pos=}')

        if not self.reference:
            raise ValueError('cannot compute DE without a reference')

        if len(self.replicates) < 2:
            print(f'DE: {self.name} has not enough replicates to compute DE.')
            return None
            # pd.DataFrame(columns=['log2FoldChange', 'log10pvalue', 'log10padj'])

        cond = self.name
        ref = self.reference.name

        df = self.reference.get_counts(pos=pos, **kwargs).join(
            self.get_counts(pos=pos, **kwargs)
        )

        if filter_size:
            df = df[~df.index.str.startswith(' ')]

        sample_df = pd.DataFrame.from_records(
            [
                {
                    'samplename': str(replicate),
                    'sample': sample.name,
                    'replicate': replicate.replicate,
                }
                for sample in [self, self.reference]
                for replicate in sample.replicates
            ]
        ).set_index('samplename', drop=False)

        _cache_prefix = f'{self.name_ref}{"_pos="+str(pos) if pos else ""}'
        if kwargs:
            _cache_prefix += f'_{"_".join([f"{key}={value}" for key, value in kwargs.items()])}'

        return DE(
            df,
            sample_df,
            cond=cond,
            ref=ref,
            join=join,
            quiet=quiet,
            multi=multi,
            n_cpus=n_cpus,
            raw=raw,
            _nocache=_nocache,
            _cache_dir=self._cache_dir,
            _cache_prefix=_cache_prefix,
        )

    def hmap(
        self,
        r=None,
        c=None,
        *,
        pos=None,
        #how='aax',
        col='auto',
        transform=np.log2,
        cmap='vlag',
        vmax=None,
        center=None,
        ax=None,
        heatmap_kwargs=None,
        **kwargs,
    ):
        """
        Generates a heatmap of enrichment for combinations of 2 positions.

        Parameters
        ----------
        r : str
            The row position on the ribosome for the heatmap.
        c : str
            The column position on the ribosome for the heatmap.
        pos : str or list
            Either a specific position in the form "r:c" or a list of positions to analyze.
        how : str
            Defines the method to compute the counts (e.g., 'mean', 'sum', 'count').
            If 'aax' is provided, sequences with stop codons in the peptide are excluded.
        col : str
            The dataset column used for computations.
        transform : callable, optional
            A function or callable to apply to the dataset before generating the heatmap.
        cmap : str or matplotlib.colors.Colormap
            The colormap to use for the heatmap visualization.
        vmax : float, optional
            The maximum value for color scaling in the heatmap.
        center : float, optional
            The midpoint value for centering the colormap.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. If not provided, a new figure and axes are created.
        heatmap_kwargs: dict
            Parameters passed to the sns.heatmap method
        kwargs : dict
            Additional parameters used to filter the dataset. This allows for fine-tuning
            of the data before generating the heatmap.

        Returns
        -------
        matplotlib.axes.Axes
            The heatmap axes object containing the visualization.
        """

        if not r and not c:
            if not pos:
                raise ValueError(
                    'Must provide "pos" or both "c" and "r" parameters.'
                )
            r, c = pos.replace(',', ':').split(':')
        pos = f'{r},{c}'
        table = self.get_counts_ratio(pos=pos, **kwargs)

        if col == 'auto':
            col = 'ratio'
            center = 0
            cmap = 'vlag'
            transform = np.log2

        how = kwargs.get('how', 'aax')
        if how in ('aa', 'aax'):
            names = r, c
            colors = aa_colors
            order = aa_order
        else:
            colors = {
                'A': '#008000',
                'C': '#0000FF',
                'G': '#FFA806',
                'T': '#FF0000',
            }
            order = list('ACGT')
            names = r, c

        table = table.copy()  # avoid modifying the input object
        # keep only motifs without gap
        try:
            table = table.loc[~table.index.str.contains(' ')]
        except:
            print(table)

        names = list(map(get_ribosome_site, names))

        table.index = pd.MultiIndex.from_tuples(
            [tuple(x) for x in table.index], names=names
        )

        if not ax:
            f, ax = plt.subplots()

        t = table[col].unstack()
        if len(t) == 0:
            return ax

        d = t.reindex(index=order, columns=order)

        if transform:
            d = transform(d)

        if not vmax:
            # vmax = max(-d.values.min(), d.values.max())
            vmax = np.nanmax(abs(d.values))

        if not heatmap_kwargs:
            heatmap_kwargs = {}
        heatmap_kwargs.setdefault('vmax', vmax)
        heatmap_kwargs.setdefault('vmin', -vmax if col == 'ratio' else 0)
        heatmap_kwargs.setdefault('xticklabels', True)
        heatmap_kwargs.setdefault('yticklabels', True)
        heatmap_kwargs.setdefault('center', center)
        heatmap_kwargs.setdefault('cmap', cmap)
        heatmap_kwargs.setdefault('square', True)
        heatmap_kwargs.setdefault('linecolor', '#CECECE')
        heatmap_kwargs.setdefault('linewidths', 0.05)
        heatmap_kwargs.setdefault('cbar_kws', {'shrink': 0.5})

        sns.heatmap(d, ax=ax, **heatmap_kwargs)

        for a in ax.figure.axes:
            if a.get_label() == '<colorbar>':
                a.set_ylabel('log2(enrichment)' if col == 'ratio' else col)
        for s in ax.spines.values():
            s.set_visible(True)

        # horizontal centered yticklabels
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=0,
            ha='center',
            fontdict={'fontfamily': 'monospace'},
        )
        ax.set_yticklabels(
            ax.get_yticklabels(),
            rotation=90,
            ha='center',
            fontdict={'fontfamily': 'monospace'},
        )

        ## AAs colors
        for a in [ax.xaxis, ax.yaxis]:
            for t in a.get_ticklabels():
                lbl = t.get_text()
                bbox = dict(
                    boxstyle='square',
                    pad=0.15,
                    ec='none',
                    fc=colors[lbl],
                    alpha=0.5,
                )
                t.set_bbox(bbox)
                t.get_bbox_patch().set_width(5)

        ax.figure.canvas.draw()

        for t in ax.get_xticklabels():
            t.get_bbox_patch().set_width(5)
        for t in ax.get_yticklabels():
            t.get_bbox_patch().set_width(5)

        ax.invert_yaxis()

        for t in ax.get_xticklabels():
            t.get_bbox_patch().set_width(5)
        for t in ax.get_yticklabels():
            t.get_bbox_patch().set_width(5)

        return ax

    def hmap_grid(
        self,
        pos=None,
        #how='aax',
        col='auto',
        transform=np.log2,
        cmap='vlag',
        vmax=None,
        center=None,
        **kwargs,
    ):
        """
        Creates a grid of heatmaps for all combinations of ribosome positions passed in `pos`.

        Each cell in the upper triangle of the grid represents a heatmap of enrichment between two positions,
        with the visualization parameters inherited from the `hmap` method.

        Parameters
        ----------
        pos : iterable, optional
            An iterable of ribosome positions for generating combinations (e.g., ['-2', 'E', 'P', 'A']).
            If not provided, defaults to the set of positions ['-2', 'E', 'P', 'A'].
        how : str, optional
            If 'aax' is provided, sequences with stop codons in the peptide are excluded.
        col : str, optional
            The dataset column used for computations. Displays the enrichment by default.
        transform : callable, optional
            A function or callable to apply to the dataset before generating the heatmaps.
            Defaults to `numpy.log2`.
        cmap : str or matplotlib.colors.Colormap, optional
            The colormap to use for the heatmap visualizations. Defaults to 'vlag'.
        vmax : float, optional
            The maximum value for color scaling in the heatmaps.
        center : float, optional
            The midpoint value for centering the colormap.
        kwargs : key, value pairings
            Additional parameters used to filter the dataset or control heatmap generation via the `hmap` method.

        Returns
        -------
        matplotlib.figure.Figure
            The figure object containing the grid of heatmaps.
        """

        if not pos:
            pos = ['-2', 'E', 'P', 'A']

        N = len(pos) - 1
        f, axes = plt.subplots(nrows=N, ncols=N, figsize=(N * 5, N * 5))

        for i, r in enumerate(pos[:-1]):
            for j, c in enumerate(pos[1:]):
                ax = axes[i, j]
                if i > j:
                    ax.set_visible(False)
                    continue
                # if i == j:
                #    ax.set_ylabel(r)
                if i == 0:
                    ax.set_title(c)
                self.hmap(
                    r=r,
                    c=c,
                    #how=how,
                    col=col,
                    transform=transform,
                    cmap=cmap,
                    vmax=vmax,
                    center=center,
                    ax=ax,
                    **kwargs,
                )
        return f

    @lru_cache
    def get_counts_ratio_pos(self, pos=None, **kwargs):
        """
        Computes a DataFrame with the enrichment ratios for each ribosome position.

        This method calculates the enrichment for amino acids at the specified positions
        on the ribosome and organizes the results into a DataFrame. Each row of the DataFrame
        corresponds to a ribosome position.

        Parameters
        ----------
        pos : iterable, optional
            An iterable of ribosome positions for which to compute enrichment ratios (e.g., ('-2', 'E', 'P', 'A')).
            If not provided, defaults to ('-2', 'E', 'P', 'A').
        how : str, optional
            If 'aax' is provided, sequences with stop codons in the peptide are excluded.
        **kwargs : dict, optional
            Additional parameters to filter the data or customize the ratio computations.

        Returns
        -------
        pandas.DataFrame
            A DataFrame where rows correspond to ribosome positions and columns correspond to amino acids
            (ordered by a predefined amino acid sequence). The values in the DataFrame represent
            the enrichment ratios for each position and amino acid.
        """

        if not pos:
            pos = ('-2', 'E', 'P', 'A')

        return (
            pd.concat(
                {
                    p: self.get_counts_ratio(p, **kwargs)['ratio']
                    for p in pos
                },
                axis=1,
            )
            .reindex(aa_order)
            .rename_axis(index='amino-acid', columns='site')
            .T
        )

    def hmap_pos(
        self,
        pos=None,
        *,
        #how='aax',
        #col='auto',
        #transform=np.log2,
        cmap='vlag',
        vmax=None,
        center=0,
        ax=None,
        **kwargs,
    ):
        """
        Generates a heatmap of enrichment ratios for amino acid positions across ribosome sites.

        This method visualizes the enrichment ratios as a heatmap, where the rows correspond
        to different ribosome positions and the columns represent amino acids.

        Parameters
        ----------
        pos : tuple, optional
            Ribosome positions for which to compute and visualize enrichment ratios (e.g., `('-2', 'E', 'P', 'A')`).
        how : str, optional
            If 'aax' is provided, sequences with stop codons in the peptide are excluded. Default is 'aax'.
        col : str, optional
            The DataFrame column to utilize for enrichment visualization. Defaults to 'auto'.
        transform : callable, optional
            A function or callable to apply to the enrichment matrix before plotting. Defaults to `numpy.log2`.
        cmap : str or matplotlib.colors.Colormap, optional
            The colormap to use for the heatmap visualization. Defaults to 'vlag'.
        vmax : float, optional
            The maximum value for color scaling in the heatmap. If not provided, it defaults to the
            maximum absolute value in the enrichment matrix.
        center : float, optional
            The midpoint of the colormap. Defaults to 0.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. A new figure and axes are created if not provided.
        **kwargs : dict, optional
            Additional parameters to customize the enrichment computation or filtering.

        Returns
        -------
        matplotlib.axes.Axes
            The axes object containing the heatmap visualization.

        Notes
        -----
        - The rows of the heatmap correspond to ribosome positions, while the columns represent amino acids.
        - Tick labels are styled using the `aa_colors` dictionary to match the biochemical categories of amino acids.
        - Enrichment ratios are automatically log2-transformed by default.
        """

        df_single = self.get_counts_ratio_pos(pos=pos, **kwargs)

        if not vmax:
            vmax = np.nanmax(df_single.abs())

        if not ax:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        sns.heatmap(
            df_single,
            yticklabels=True,
            xticklabels=True,
            center=center,
            cmap=cmap,
            square=True,
            linecolor='#CECECE',
            linewidths=0.05,
            cbar_kws={
                'shrink': 0.75,
                # 'ticks': [-0.5,-0.25,0,0.25,0.5],
                # 'ticks': [-1,-0.5,0,0.5,1],
                'orientation': 'horizontal',
                'label': 'log2(enrichment)',
            },
            vmax=vmax,
            vmin=-vmax,
            ax=ax,
        )
        for s in ax.spines.values():
            s.set_visible(True)

        for a in [ax.xaxis]:
            for t in a.get_ticklabels():
                aa = t.get_text()
                bbox = dict(
                    boxstyle='square',
                    pad=0.15,
                    ec='none',
                    fc=aa_colors[aa],
                    alpha=0.5,
                )
                t.set_bbox(bbox)
                t.get_bbox_patch().set_width(5)

        ax.set_yticklabels(
            ax.get_yticklabels(),
            rotation=0,
            ha='center',
            fontdict={'fontfamily': 'monospace'},
        )
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=0,
            ha='center',
            fontdict={'fontfamily': 'monospace'},
        )
        fig.canvas.draw()
        [t.get_bbox_patch().set_width(5) for t in ax.get_xticklabels()]

        # fig.tight_layout()
        return ax

    def volcano(
        self,
        pos=None,
        query='auto',
        motif=None,
        ax=None,
        x='log2FoldChange',
        y='log10pvalue',
        query_color='#BC0909',
        motif_color='#EDEA20',
        color='k',
        params={'s': 2, 'figsize': (6, 4.5)},
        annotate=True,
        text_stroke=False,
        outfile=None,
        density_thresh=0,
        **kwargs,
    ):

        df = self.DE(pos=pos, **kwargs)
        if df is None:
            print(f'Skipping volcano plot for {self}')
            return
        # df.astype(float)  ### FIXME

        df['rank'] = df[x].abs() * df[y]
        df['qrank'] = df['rank'].rank(pct=True)
        df.sort_values(by='qrank', ascending=False, inplace=True)

        if not ax:
            f, ax = plt.subplots()

        if not annotate:
            annotate = ''

        df.plot.scatter(x=x, y=y, color=color, ax=ax, **params)

        dfs = []
        colors = []

        if m := re.search(r'head:\s*(\d+)', query):
            head = df['qrank'].head(int(m.group(1))).iloc[-1]
            query = f'qrank > {head}'

        elif m := re.search(r'density:\s*(\d+)', query):
            df = df.assign(
                keep=flag_low_density(df[x], df[y], 50, int(m.group(1)))
            )
            query = f'keep'

        elif query.startswith('auto'):
            quantile = 0.90
            m = re.search(r'auto:\s*(\d+\.?\d*)', query)
            if m:
                quantile = float(m.group(1))
                query = 'auto'
            if query == 'auto':
                x_thresh = df['log2FoldChange'].abs().quantile(quantile)
                y_thresh = df['log10pvalue'].abs().quantile(quantile)

                x_min = df['log2FoldChange'].agg(['min', 'max']).abs().min()
                x_thresh = min(x_thresh, x_min)

                query = f'(log10pvalue > @y_thresh) & (abs(log2FoldChange) > @x_thresh)'
            else:
                thresh = {2: (0.05, 0.2), 3: (0.05, 1), 4: (0.2, 2)}
                padj, l2fc = thresh.get(len(df.index[0]), (0.05, 1))
                query = f'(padj < {padj}) & (log2FoldChange > {l2fc})'

        if query:
            df_query = df.query(
                query, engine='python'
            )  # engine='python' to prevent bug with <NA> type
            df_query.plot.scatter(x=x, y=y, color=query_color, ax=ax, **params)
            if (annotate == True) or ('query' in annotate):
                dfs.append(df_query)
                colors.append(query_color)

        if motif:
            df_motif = df[df.index.str.match(motif)]
            df_motif.plot.scatter(x=x, y=y, color=motif_color, ax=ax, **params)
            if (annotate == True) or ('motif' in annotate):
                dfs.append(df_motif)
                colors.append(motif_color)

        if annotate:
            from matplotlib import patheffects

            for d, c in zip(dfs, colors):
                if density_thresh:
                    d = d[
                        flag_low_density(
                            d[x], d[y], bins=30, thresh=density_thresh
                        )
                    ]
                for name, sx, sy in zip(
                    d.index, d[x], d[y]
                ):   # should be much faster than iterrows if many points

                    # adapt position of the label
                    left = ((sx < 0) and (sy < -sx)) or (
                        (sx > 0) and (sy > sx)
                    )

                    txt = ax.annotate(
                        name,
                        (sx, sy),
                        xytext=(-5 if left else 5, 0),
                        textcoords='offset pixels',
                        ha='right' if left else 'left',
                        va='center',
                    )
                    if text_stroke:
                        stroke_color = (
                            c if isinstance(text_stroke, bool) else text_stroke
                        )
                        txt.set_path_effects(
                            [
                                patheffects.withStroke(
                                    linewidth=2, foreground=stroke_color
                                )
                            ]
                        )

        # make sure the Y-axis is at least showing the [0, 1] range for -log10(p-value)
        # this avoid having a weird graph is nothing is significant
        if y.startswith('log10p'):
            ax.set_ylim(top=max(1.1, ax.get_ylim()[1]))

        ax.set_title(f'{self.name_vs_ref}: {pos}')

        if outfile:
            ax.figure.savefig(outfile)

        return ax

    def all_logos(self, logo_kwargs=None, **kwargs):
        f, axes = plt.subplots(nrows=len(self.replicates), sharex=True)
        for r, ax in zip(self.replicates, axes.flat):
            r.logo(logo_kwargs=logo_kwargs, ax=ax, **kwargs)
            ax.set_title(r)

    def logo(self, pos=None, logo_kwargs=None, ax=None, **kwargs):

        df = np.log2(
            self.get_counts_ratio_pos(pos=pos, **kwargs)
        ).fillna(0)

        import logomaker

        if logo_kwargs is None:
            logo_kwargs = {'color_scheme': 'NajafabadiEtAl2017'}

        # df = logomaker.transform_matrix(df, from_type='counts', to_type=type)
        # print(df)

        if not ax:
            ax = plt.subplot()

        logo = logomaker.Logo(df.reset_index(drop=True), ax=ax, **logo_kwargs)
        logo.ax.set_xticks(range(df.shape[0]), df.index)
        return logo

    @cached_property
    def itp_len(self):
        """
        Combines the counts of inverse-toeprints (ITPs) for each length across all replicates.

        This method extracts the counts of inverse-toeprints for each length from the metadata of each
        replicate and combines them into a single DataFrame, keeping the data for each replicate independent.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with the following columns:
            - `length` : int
                The length of the inverse-toeprints.
            - `replicate` : str
                The replicate identifier.
            - `count` : int
                The count of inverse-toeprints of the given length for the replicate.
            - `sample` : str
                The name of the sample this data belongs to.
        """

        return (
            pd.concat(
                {
                    r.replicate: pd.Series(r.meta['fastq']['lengths'])
                    for r in self.replicates
                },
                axis=1,
            )
            .rename(int)
            .rename_axis('length')
            .reset_index()
            .melt('length', var_name='replicate', value_name='count')
            .assign(sample=self.name)
        )

    def itp_len_plot(
        self, ax=None, min_codon=0, max_codon=10, limit=100, norm=False
    ):
        """
        Generates a line plot of inverse-toeprint (ITP) counts per length.

        This method uses the output of `itp_len` to create a line plot showing the counts of
        inverse-toeprints across lengths for each replicate. Optionally, counts can be normalized
        (per million reads), and the plotted lengths can be limited.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes to draw the plot on. A new figure and axes are created if not provided.
        min_codon : int, optional
            The minimum codon position to annotate on the plot. Defaults to 0.
        max_codon : int, optional
            The maximum codon position to annotate on the plot. Defaults to 10.
        limit : int, optional
            The maximum length to include in the plot. Defaults to 100.
        norm : bool, optional
            Whether to normalize counts to reads per million. Defaults to False.

        Returns
        -------
        matplotlib.axes.Axes
            The axes object containing the plotted lineplot.

        Notes
        -----
        - The x-axis represents the distance from the 3' end of the inverse-toeprint in nucleotides.
        - The y-axis shows the counts of inverse-toeprints, either absolute or normalized per million reads.
        - Each replicate is plotted independently and distinguished by the `hue` attribute in the plot.
        """

        df_len = self.itp_len.copy()

        names = {
            'count': 'Number of reads (per million)'
            if norm
            else 'Number of reads',
            'length': 'Distance from the start site [nt]',
            'rep': 'Replicate',
            'sample': 'sample',
        }

        if limit:
            df_len = df_len.loc[df_len['length'].le(limit)]

        if norm:
            df_len['count'] /= (
                df_len.groupby(['sample', 'replicate'])['count']
                .transform('sum')
                .div(1_000_000)
            )

        if not ax:
            f, ax = plt.subplots()

        sns.lineplot(
            data=df_len.rename(columns=names),
            x=names['length'],
            y=names['count'],
            hue='replicate',
            ax=ax,
        )
        itp_len_add_positions(ax, min_codon=min_codon, max_codon=max_codon)

        return ax

    @cached_property
    def toeprint_df(self):
        return (
            pd.concat(
                {
                    r.name: pd.Series(r.meta['fastq']['lengths'])
                    for r in self.replicates
                },
                axis=1,
            )
            .rename(int)
            .sort_index(ascending=False)
        )

    def itoeprint(
        self,
        plot='bands',
        norm='mean',
        norm_range=(21, 51),
        exposure=1,
        limit=(0, 100),
        show_range=False,
        ax=None,
        interactive=False,
    ):
        if interactive:
            from ipywidgets import interact, fixed
            import ipywidgets as widgets

            interact(
                itoeprint_plot,
                dataset=fixed(self),
                plot=widgets.Dropdown(options=['bands', 'shades'], value=plot),
                norm=widgets.Dropdown(
                    options=[None, 'mean', 'median', 'max', 'std'], value=norm
                ),
                norm_range=widgets.IntRangeSlider(
                    value=norm_range, min=0, max=100
                ),
                exposure=widgets.FloatSlider(
                    value=exposure, min=0.1, max=10.0, step=0.1
                ),
                limit=widgets.IntRangeSlider(value=limit, min=0, max=500),
                show_range=widgets.Checkbox(value=True),
                # equalize_lanes=widgets.Checkbox(value=equalize_lanes),
                ax=fixed(None),
            )
        else:
            return itoeprint_plot(
                self,
                plot=plot,
                norm=norm,
                norm_range=norm_range,
                exposure=exposure,
                limit=limit,
                show_range=show_range,
                ax=ax,
            )


class DataSet:
    r"""
    Loads an iTP-Seq dataset and provides methods for analyzing and visualizing the data.

    A DataSet object is constructed to handle iTP-Seq Samples with their respective Replicates.
    By default, it infers the files to uses in the provided directory by looking for "*.processed.json" files produced
    during the initial step of pre-processing and filtering the fastq files.
    It uses the pattern of the file names to group the Replicates into a Sample, and to define which condition
    is the reference in the DataSet (the Sample with name "noa" by default).

    Attributes
    ----------
    data_path : str or Path
        Path to the data directory containing the output files from the fastq pre-processing.
    result_path: str or Path
        Path to the directory where the results of the analysis will be saved.
    samples: dict or None
        Dictionary of Samples in the DataSet. By default, it is None and will be populated automatically.
    keys: tuple
        Properties in the file name to use for identifying the reference.
    ref_labels: str or tuple
        Specifies the reference: e.g. 'noa' or (('sample', 'noa'),)
    cache_path: str or Path
        Path used to cache intermediate results. By default, this creates a subdirectory called "cache" in the result_path directory.
    file_pattern: str
        Regex pattern used to identify the sample files in the data_path directory.
        If None, defaults to `r'(?P<lib_type>[^_]+)_(?P<sample>[^_\d]+)(?P<replicate>\d+)\.processed\.json'`
        which matches files like nnn15_noa1.processed.json, nnn15_tcx2.processed.json, etc.
    aafile_pattern: str
        Pattern used to identify the amino acid files in the data_path directory.
        It will use the values captured in the `file_pattern` regex to construct the file names.
        If None, defaults to `'{lib_type}_{sample}{replicate}_aa.processed.txt'`

    Examples
    --------
    Creating a DataSet from a simple antibiotic treatment (tcx) vs no treatement (noa) with 3 replicates each (1, 2, 3).

    Load a dataset from the current directory, inferring the samples automatically
     >>> from itpseq import DataSet
     >>> data = DataSet(data_path='.')
     >>> data
     DataSet(data_path=PosixPath('.'),
        reference=Sample(noa:[1, 2, 3]),
        samples=[Sample(noa:[1, 2, 3]),
                 Sample(tcx:[1, 2, 3], ref: noa)],
        )

    Compute a standard report and export it as PDF
     >>> data.report('my_experiment.pdf')

    Display a graph of the inverse-toeprints lengths for each sample
     >>> data.itp_len_plot(row='sample')
    """

    def __init__(
        self,
        data_path: Path = '.',
        result_path: Path = None,
        samples: Optional[dict] = None,
        keys=None,
        ref_labels: Optional[Union[str, tuple]] = 'noa',
        # max_workers: int = 4,
        cache_path=None,
        file_pattern=None,
        aafile_pattern=None,
    ):
        self.data_path = Path(data_path)
        if result_path is None:
            result_path = self.data_path / 'results'
        if cache_path is None:
            cache_path = self.data_path / 'results/cache'
        self.result_path = Path(result_path)
        if not self.result_path.exists():
            self.result_path.mkdir(parents=True)
        self.cache_path = Path(cache_path)
        if not self.cache_path.exists():
            self.cache_path.mkdir(parents=True)
        # self.file_pattern = r'nnn15_(?P<sample>[^_]+)(?P<replicate>\d+)\.processed\.json'
        self.file_pattern = (
            file_pattern
            or r'(?P<lib_type>[^_]+)_(?P<sample>[^_\d]+)(?P<replicate>\d+)\.processed\.json'
        )
        self.aafile_pattern = (
            aafile_pattern or '{lib_type}_{sample}{replicate}_aa.processed.txt'
        )
        # if no keys are provided, use all the named capturing groups of file_pattern, excepted "replicate"
        self.keys = (
            keys
            if keys
            else [
                g
                for g in re.compile(self.file_pattern).groupindex
                if g != 'replicate'
            ]
        )

        # ref_labels can be passed as various formats
        # 'noa'
        # (('sample', 'noa'),)
        # (('sample', 'noa'), ('lib_type', 'nnn15'))
        # {'sample': 'noa'}
        # and will be converted to a tuple of (key, value) tuples
        if isinstance(ref_labels, str):
            self.ref_labels = (('sample', ref_labels),)
        elif isinstance(ref_labels, dict):
            self.ref_labels = dict_to_tuple(ref_labels)
        elif isinstance(ref_labels, tuple):
            self.ref_labels = ref_labels
            assert (
                isinstance(ref_labels, tuple)
                and (len(ref_labels) == 0)
                or all(isinstance(x, tuple) for x in ref_labels)
            )
        else:
            self.ref_labels = ()

        self.reference = None
        # self.max_workers = max_workers

        self.samples = samples or self._infer_samples()
        self.replicates = {
            r.name: r for s in self.samples.values() for r in s.replicates
        }

    def __getitem__(self, key):
        try:
            return self.samples[key]
        except KeyError:
            return self.replicates[key]

    def __repr__(self):
        sep = ',\n        '
        sep2 = sep + ' ' * 9
        filepat = f'{sep}file_pattern={self.file_pattern!r}'
        samples = (
            f'{sep}samples=[{sep2.join(repr(s) for s in self.samples.values())}]{sep}'
            if self.samples
            else ''
        )
        ref = (
            f'{sep}reference={repr(self.reference)}' if self.reference else ''
        )
        # return f'''DataSet(data_path={repr(self.data_path)}, result_path={repr(self.result_path)}{ref}{samples})'''
        return f"""DataSet(data_path={repr(self.data_path)}{filepat}{ref}{samples})"""

    @property
    def samples_with_ref(self):
        return {k: s for k, s in self.samples.items() if s.reference}

    def _infer_samples(self) -> list[str]:
        """Infers sample names from the files in the data path."""
        inferred_samples = defaultdict(list)
        # file_paths = list(self.data_path.glob("*.json"))

        for f in self.data_path.iterdir():  ## TODO: wrap with SampleGrouper
            if m := re.search(self.file_pattern, f.name):
                labels = m.groupdict()
                fastq_meta = read_log_json(f)
                # if set(labels) >= {'sample', 'replicate'}:
                if set(labels) > {'replicate'}:
                    rep = labels['replicate']
                    # keys = tuple(sorted(t for t in labels.items() if t[0] in self.keys))
                    keys = dict_to_tuple(labels, keep=self.keys)
                    inferred_samples[keys].append(
                        {
                            'filename': self.data_path
                            / self.aafile_pattern.format(**labels),
                            'replicate': rep,
                            'labels': labels,
                            'fastq': fastq_meta,
                        }
                    )

        samples = {}
        samples_by_tup = {}
        # self.infsamp = inferred_samples
        # print(inferred_samples)

        for key, data in inferred_samples.items():
            s = Sample(
                labels=dict(key), data=data, dataset=self, keys=self.keys
            )
            samples[s.name] = s
            samples_by_tup[key] = s

        for s in samples.values():
            # get potential ref key/label
            ref_lbl = s.labels | dict(self.ref_labels)
            ref_tup = tuple(sorted(ref_lbl.items()))
            if ref_lbl != s.labels and ref_tup in inferred_samples:
                s.reference = samples_by_tup[ref_tup]

        return {k: samples[k] for k in sorted(samples)}

    def reorder_samples(self, order, validate=True, reorder_replicates=True):
        if validate:
            assert set(order) == set(self.samples.keys())
        self.samples = {k: self.samples[k] for k in order}
        if reorder_replicates:
            keys = {k: i for i, k in enumerate(order)}
            self.replicates = {
                r.name: r
                for r in sorted(
                    self.replicates.values(),
                    key=lambda r: keys.get(r.sample_name),
                )
            }
            # force update cached attributes that depend on a specific order of the replicates
            for attr in ['toeprint_df']:
                if attr in self.__dict__:
                    self.__dict__.pop(attr)
        # force update cached attributes that depend on a specific order of the samples
        # self.__dict__.pop('toeprint_df')

    def infos(self, html=False):
        """Displays summary information about the dataset sequences.

        Atttributes
        -----------
        html: bool
            if True, returns the table as HTML, otherwise as DataFrame (default)

        """
        order = [
            'total_sequences',
            'noadaptor',
            'contaminant',
            'lowqual',
            'tooshort',
            'toolong',
            'extra0',
            'extra1',
            'extra2',
        ]
        key = {c: i for i, c in enumerate(order)}
        out = (
            pd.DataFrame.from_dict(
                {name: r.meta['fastq'] for name, r in self.replicates.items()},
                orient='index',
            )
            .sort_index(axis=1, key=lambda x: x.map(key))
            .select_dtypes('number')
        )
        if html:
            return table_to_html(out)
        return out

    def DE(self, pos='E:A', **kwargs):
        """Computes the log2-FoldChange for each motif described by `pos` for each sample in the DataSet relative to their reference

        Atttributes
        -----------
        pos: str
            position of the motif to consider.
            This ca be a range of positions (e.g. '-2:A' for -2/E/P/A sites)
            or a combination of disjoint positions (e.g. 'E,A' for the combination of E and A sites).
        **kwargs: dict
            parameters passed to `Sample.get_counts` computes the counts for the given motif.
            for example `min_peptide=3` to consider only peptides of at least 3 amino acids.

        ."""
        out = {}
        for sample in self.samples_with_ref.values():
            # out[(sample.name, pos)] = sample.DE(pos=pos, **kwargs)
            out[sample.name] = sample.DE(pos=pos, **kwargs)
        return out

    def itp_len_plot(
        self,
        ax=None,
        col=None,
        row=None,
        min_codon=0,
        max_codon=10,
        agg=False,
        limit=100,
        norm=True,
        hue='auto',
        plt_kwargs=dict(kind='line', height=2, aspect=3),
    ):
        df_len = pd.concat(
            (s.itp_len for s in self.samples.values()), ignore_index=True
        )

        names = {
            'count': 'Number of reads (per million)'
            if norm
            else 'Number of reads',
            'length': 'Distance from the start site [nt]',
            'rep': 'Replicate',
            'sample': 'sample',
        }

        if limit:
            df_len = df_len.loc[df_len['length'].le(limit)]

        if norm:
            df_len['count'] /= (
                df_len.groupby(['sample', 'replicate'])['count']
                .transform('sum')
                .div(1_000_000)
            )

        if col or row:
            g = sns.relplot(
                data=df_len.rename(columns=names),
                col=col,
                row=row,
                x=names['length'],
                y=names['count'],
                hue='replicate' if hue == 'auto' else hue,
                **plt_kwargs,
            )
            for ax in g.axes.flat:
                itp_len_add_positions(
                    ax, min_codon=min_codon, max_codon=max_codon
                )
            return g

        if not ax:
            f, ax = plt.subplots()

        sns.lineplot(
            data=df_len.rename(columns=names),
            x=names['length'],
            y=names['count'],
            hue='sample' if hue == 'auto' else hue,
            ax=ax,
        )
        itp_len_add_positions(ax, min_codon=min_codon, max_codon=max_codon)

        return ax

    @cached_property
    def toeprint_df(self):
        return (
            pd.concat(
                {
                    r.name: pd.Series(r.meta['fastq']['lengths'])
                    for _, r in self.replicates.items()
                },
                axis=1,
            )
            .rename(int)
            .sort_index(ascending=False)
        )

    itoeprint = Sample.itoeprint

    def report(self, template='report', output=None):
        from jinja2 import Environment, FileSystemLoader, PackageLoader

        template_dir = Path(__file__).resolve().parent / 'templates'
        static_dir = Path(__file__).resolve().parent / 'static'
        env = Environment(loader=FileSystemLoader(template_dir))
        if not template_dir.joinpath(f'{template}.html').exists():
            raise ValueError(
                f'Template {template}.html not found in {template_dir}'
            )
        html_template = env.get_template(f'{template}.html')

        # additional variables
        vars = {
            'today': datetime.date.today(),
        }

        if output:
            output = Path(output)
            if output.suffix == '.html':
                css_file = static_dir / f'{template}-html.css'
                if not css_file.exists():
                    css_file = static_dir / f'report-html.css'
                html = html_template.render(
                    dataset=self,
                    plot_to_html=plot_to_html,
                    css_file=css_file,
                    output='html',
                    **vars,
                )
                with open(output, 'w') as f:
                    f.write(html)
            elif output.suffix == '.pdf':
                from weasyprint import HTML, CSS

                css_file = static_dir / f'{template}-pdf.css'
                if not css_file.exists():
                    css_file = static_dir / f'report-pdf.css'
                html = html_template.render(
                    dataset=self,
                    plot_to_html=plot_to_html,
                    output='pdf',
                    **vars,
                )
                HTML(string=html).write_pdf(
                    output, stylesheets=[CSS(filename=css_file)]
                )
            return
        else:
            html = html_template.render(
                dataset=self, plot_to_html=plot_to_html, output='other', **vars
            )

            if IPYTHON:
                from IPython.display import HTML

                return HTML(html)

            return html


if __name__ == '__main__':
    pass
