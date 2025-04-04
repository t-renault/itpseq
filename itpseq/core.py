"""Core classes for itpseq"""

import datetime
import re
import uuid
import warnings
from collections import defaultdict
from copy import deepcopy
from functools import cached_property, lru_cache, total_ordering, wraps
from pathlib import Path
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas.api.extensions import no_default

from .config import *
from .plotting import *
from .processing import *
from .utils import *

IPYTHON = False
try:
    from IPython import get_ipython

    IP = get_ipython()
    if IP is not None and getattr(IP, 'kernel', None) is not None:
        IPYTHON = True
        del IP
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
        file_prefix: Optional[Path] = None,
        sample: Optional['Sample'] = None,
        labels: Optional[dict] = None,
        **kwargs,
    ):
        # print(f'Replicate {kwargs=}')
        self.name = None
        self.file_prefix = file_prefix
        self.sample = sample
        self.replicate = replicate
        self.labels = labels
        self.rename()   # automatically set replicate name
        self._cache_dir = self.dataset.cache_path if self.dataset else None
        self.meta = kwargs

    def __repr__(self):
        return f'Replicate({self.name})'

    def __str__(self):
        return self.name

    def __gt__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        if not self.sample and other.sample:
            try:
                return self.replicate > other.replicate
            except TypeError:
                return str(self.replicate) > str(other.replicate)
        else:
            try:
                return (self.sample, self.replicate) > (
                    other.sample,
                    other.replicate,
                )
            except TypeError:
                return (self.sample, str(self.replicate)) > (
                    other.sample,
                    str(other.replicate),
                )

    @property
    def filename(self):
        """File name for the nucleotide inverse toeprints for the replicate"""
        return self.file_prefix + f'.nuc.{ITP_FILE_SUFFIX}.txt'

    @property
    def aa_filename(self):
        """File name for the amino-acid inverse toeprints for the replicate"""
        return self.file_prefix + f'.aa.{ITP_FILE_SUFFIX}.txt'

    @property
    def json_filename(self):
        """File name for the JSON metadata for the replicate"""
        return self.file_prefix + f'.{ITP_FILE_SUFFIX}.json'

    @property
    def sample_name(self):
        """Linked Sample name (if any)"""
        return self.sample.name if self.sample else ''

    @property
    def dataset(self):
        """Linked DataSet name (if any)"""
        return self.sample.dataset if self.sample else None

    def rename(self, name=None):
        """
        Sets the name of the replicate from a parameter or automatically.

        Parameters
        ----------
        name : str, optional
            name to use as the new name for the replicate.

        Examples
        --------
        Rename the replicate with a parameter
         >>> rep.rename(name='new_name')

        Rename the replicate automatically from its parent sample data
         >>> rep.rename()
        """
        if name is None:
            self.name = (
                f'{self.sample.name}.{self.replicate}'
                if self.sample and self.replicate
                else self.replicate
            )
        else:
            self.name = name

    def copy(self, sample=None, replicate=None, name=None):
        """
        Creates a copy of the replicate.

        Parameters
        ----------
        sample : Sample, optional
            New sample for the replicate.
        replicate : optional
            New identifier or label for the replicate.
        name : str, optional
            Instead of passing a Sample and replicate identifier, use a string as new name for the replicate.
            This will not assign a reference to a Sample.

        Returns
        -------
        Replicate
            A new Replicate object with the same data as the original Replicate and optionally an updated name.

        Examples
        --------
        Create a copy of "replicate" and change its name to "new_name".
         >>> new_rep = replicate.copy(replicate='new_name')
         >>> print(new_rep)
         Replicate(new_name)
        """
        new = deepcopy(self)
        if sample or replicate:
            new.sample = sample
            new.replicate = replicate
            new.rename()
        else:
            new.rename(name)
        return new

    def infos(self):
        """Returns the Replicate metadata from the JSON file"""
        return read_log_json(self.json_filename)

    @wraps(read_itp_file_as_series)
    def load_data(self, how='aax', **kwargs):
        """
        Reads the inverse toeprint sequences from the appropriate file
        Delegated to processing.read_itp_file_as_series

        """
        if how in {'aa', 'aax'}:
            filename = self.aa_filename
        # elif how == 'codons':  # TODO
        #    filename = self.filename
        elif how == 'nuc':
            filename = self.filename
        else:
            raise ValueError(
                f"how should be among ['aa', 'aax', 'nuc'], received: {how}"
            )

        return read_itp_file_as_series(
            filename,
            how=how,
            _cache_dir=self._cache_dir,
            _cache_prefix=self.name,
            **kwargs,
        )

    @lru_cache
    def get_counts(self, pos=None, how='aax', **kwargs):
        """
        Counts the number of reads for each motif or combination of amino-acid/position.

        Parameters
        ----------
        pos : str, optional
            Position to consider when counting the reads.
            If None is passed, then this returns a DataFrame with the counts of each amino-acid per position.
        how : str, optional
            Type of inverse toeprints to analyze (see :meth:`Replicate.load_data`).

        kwargs : optional
            Optional parameters to pass to :meth:`Replicate.load_data` (min_peptide, max_peptide, how, limit, sample)

        Returns
        -------
        Series or DataFrame
            Returns a DataFrame is pos is None, otherwise a Series.

        Examples
        --------
        Count the number of reads for each amino-acid/position combination
         >>> replicate.get_counts()
                    -8         -7         -6  ...        -1         0         1
             2879961.0  2658485.0  2449526.0  ...  793143.0   52640.0       NaN
         *         NaN        NaN        NaN  ...       NaN       NaN  910137.0
         A         NaN    12240.0    25225.0  ...  111369.0  134995.0  107591.0
         ..        ...        ...        ...  ...       ...       ...       ...
         W         NaN     2686.0     5059.0  ...   17643.0   28095.0   21577.0
         Y         NaN     9522.0    19296.0  ...   69671.0   81462.0   93099.0
         m    197624.0   221476.0   208959.0  ...  409289.0  740503.0   52640.0
         [23 rows x 10 columns]

        Count the number of reads for each motif in the E-P-A sites
         >>> replicate.get_counts(pos='E:A')
        """
        prefix = (
            '_' + '_'.join(f'{k}={v}' for k, v in sorted(kwargs.items()))
            if kwargs
            else ''
        )
        return compute_counts(
            self.load_data(**kwargs, how=how),
            pos=pos,
            how=how,
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
        Simple logo plot with default settings
         >>> logo = obj.logo()

        Logo plot with min_peptide filtering
         >>> logo = obj.logo(min_peptide=3)

        Logo plot with custom transformation type and filtering
         >>> logo = obj.logo(type='probability', min_peptides=2, fMet=True)
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

    Examples
    --------
    Get a Sample from a DataSet
     >>> sample = dataset['sample_name']

    Compute the differential expression for positions E-P-A.
     >>> sample.DE('E:A')
    """

    def __init__(
        self,
        replicates=None,
        *,
        labels=None,
        reference=None,
        dataset=None,
        # data=None,
        keys=('sample',),
        name=None,
        **kwargs,
    ):
        """"""

        self.labels = labels
        if name:
            self.name = name
        elif labels:
            # self.name = '.'.join(labels[k] or '' for k in keys) # FIXME decide of best way to handle empty keys
            self.name = '.'.join(
                labels[k] for k in keys if labels[k] is not None
            )
        else:
            self.name = f'{uuid.uuid4()}'  # if no name was provided use a UUID
        self.dataset = dataset
        self._cache_dir = self.dataset.cache_path if self.dataset else None
        self.reference = reference  # If None, this sample has no reference
        self.meta = kwargs
        # self.data = data

        if replicates is None:
            replicates = {}
        elif isinstance(replicates, list):
            names = [
                x.get('replicate', f'rep{i}')
                if isinstance(x, dict)
                else x.replicate or f'rep{i}'
                if isinstance(x, Replicate)
                else f'rep{i}'
                for i, x in enumerate(replicates, start=1)
            ]
            replicates = dict(zip(dedup_names(names), replicates))

        self.replicates = [
            Replicate(**{**data, 'replicate': k}, sample=self)
            if isinstance(data, dict)
            else data.copy(sample=self, replicate=k)
            for k, data in replicates.items()
        ]

    def __getitem__(self, key):
        return self.replicates[key]

    def __repr__(self):
        r = (
            f':[{", ".join([str(r.replicate) for r in self.replicates])}]'
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
        if other.reference == self:
            return False
        return (self.dataset, self.name) > (other.dataset, other.name)

    @property
    def name_vs_ref(self):
        """Name of the Sample combined with its reference"""
        return (
            f'{self.name} vs {self.reference.name}'
            if self.reference
            else self.name
        )

    @property
    def name_ref(self):
        """Name of the Sample combined with its reference"""
        return (
            f'{self.name}-{self.reference.name}'
            if self.reference
            else self.name
        )

    def rename(self, name, rename_replicates=True):
        r"""
        Changes the name of the sample.

        Parameters
        ----------
        name : str
            name to use as the new name for the sample.
        rename_replicates : bool
            If True (default), also rename the replicates based on the new sample name.

        Examples
        --------
        Rename the sample to "new_name".
         >>> print(sample)
         Sample(sample:[1, 2, 3], ref: noa)
         >>> sample.rename(name='new_name')
         >>> print(sample)
         Sample(new_name:[1, 2, 3], ref: noa)
         >>> print(sample.replicates)
         [Replicate(new_name.1), Replicate(new_name.2), Replicate(new_name.3)]
        """
        self.name = name
        if rename_replicates:
            for replicate in self.replicates:
                replicate.rename()

    def copy(self, name=None, reference=no_default):
        """
        Creates a copy of the sample.

        Parameters
        ----------
        name : str, optional
            New name for the sample.
        reference : Sample or None, optional
            If a parameter is used, this will set it as the reference sample.

        Returns
        -------
        Sample
            A new Sample object with the same data as the original sample and optionally an updated name and reference.

        Examples
        --------
        Create a copy of "sample" and change its name to "new_name".
         >>> new_sample = sample.copy(name='new_name')
         >>> print(new_sample)
         Sample(new_name:[1, 2, 3], ref: ref_name)

        Create a copy of "sample" with "sample2" as reference.
         >>> new_sample = sample.copy(reference=sample2)
         >>> print(new_sample)
         Sample(new_name:[1, 2, 3], ref: sample2)
        """
        new = deepcopy(self)
        if name:
            new.rename(name)
        if reference is not no_default:
            if isinstance(reference, (Sample, type(None))):
                new.reference = reference
            else:
                warnings.warn(
                    f'reference must be a Sample or None, received: {type(reference)}. '
                    f'Keeping the original reference {self.reference.name if self.reference else None}.'
                )
        return new

    def load_replicates(self, how='aax'):   ## FIXME is this useful to keep?
        """
        Loads all the Replicates in the Sample with the defined method (``how``)

        This function can be used to generate the cached files
        """
        for replicate in self.replicates:
            replicate.load_data(how=how)

    def infos(self, html=False):
        """
        Returns a table with information on the NGS reads per replicate.

        Parameters
        ----------
        html : bool
            if True, returns the table as HTML, otherwise as DataFrame (default).

        Examples
        --------
         >>> sample.infos()
                   total_sequences  noadaptor  contaminant  lowqual  tooshort  toolong   extra0   extra1   extra2  MAX_LEN
         sample.1          8384889     714414         1192   685017    385537  3341987  2308291  2528638  2833546       44
         sample.2          9120203     498202         1659   513308    104071  5664107  2673062  2972850  2976089       44
         sample.3          8490958    1043590         1328   409697    187746  4004073  2243720  2555783  2647865       44
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
                {r.name: r.infos() for r in self.replicates},
                orient='index',
            )
            .sort_index(axis=1, key=lambda x: x.map(key))
            .select_dtypes('number')
        )
        if html:
            return table_to_html(out)
        return out

    @lru_cache
    def get_counts(self, pos=None, how='aax', **kwargs):
        """
        Counts the number of reads for each motif or combination of amino-acid/position for each replicate in the sample.

        Parameters
        ----------
        pos : str, optional
            Position to consider when counting the reads.
            If None is passed, then this returns a DataFrame with the counts of each amino-acid per position.
        how : str, optional
            Type of inverse toeprints to analyze (see :meth:`Replicate.load_data`).
        kwargs : optional
            Optional parameters to pass to :meth:`Replicate.load_data` (min_peptide, max_peptide, how, limit, sample).

        Returns
        -------
        DataFrame
            Returns a DataFrame. If pos is None the columns will be a MultiIndex.

        Examples
        --------
        Count the number of reads for each amino-acid/position combination
         >>> sample.get_counts()
              sample.1                        ...  sample.3
                    -8         -7         -6  ...        -1         0         1
             2879961.0  2658485.0  2449526.0  ...  724998.0   34748.0       NaN
         *         NaN        NaN        NaN  ...       NaN       NaN  880568.0
         A         NaN    12240.0    25225.0  ...   92225.0  115164.0   85132.0
         ..        ...        ...        ...  ...       ...       ...       ...
         W         NaN     2686.0     5059.0  ...   14313.0   23730.0   17656.0
         Y         NaN     9522.0    19296.0  ...   57431.0   69162.0   81430.0
         m    197624.0   221476.0   208959.0  ...  375644.0  690250.0   34748.0
         [23 rows x 30 columns]

        Count the number of reads for each motif in the E-P-A sites
         >>> sample.get_counts(pos='E:A')
              sample.1  sample.2  sample.3
          m*  254850.0  107060.0  258338.0
          mS   54993.0   20419.0   50959.0
           m   52640.0   17860.0   34748.0
          ..        ...       ...       ...
         WFW       NaN       2.0       NaN
         WWW       NaN       1.0       NaN
         MMW       NaN       NaN       1.0
         [8842 rows x 3 columns]
        """
        # print(locals())
        return pd.concat(
            {
                str(r): r.get_counts(pos, how=how, **kwargs)
                for r in self.replicates
            },
            axis=1,
        )

    # @lru_cache
    def get_counts_ratio(
        self,
        pos=None,
        factor=1_000_000,
        exclude_empty=True,
        how='aax',
        **kwargs,
    ):
        """
        Outputs the result of `get_counts` for the sample and its reference and add extra columns:
        the normalized averages and the sample/reference ratio.

        The average is normalized for a fixed number of counts (1 million by default).

        If the sample does not have a reference, this will not compute a ratio.

        Parameters
        ----------
        pos : str, optional
            Position to consider when counting the reads (see `get_counts`).
        factor : float, optional
            The number of reads used to normalize the counts.
        exclude_empty : bool, optional
            Exclude the rows with incomplete peptides.
        how : str, optional
            Type of inverse toeprints to analyze (see :meth:`Replicate.load_data`).
        kwargs : optional
            Optional parameters to pass to :meth:`Replicate.load_data` (min_peptide, max_peptide, limit, sample).

        Returns
        -------
        DataFrame

        Examples
        --------
        Get the counts, average counts and ratio for each motif in the E-P-A sites
         >>> sample.get_counts_ratio(pos='E:A')
                 noa.1     noa.2     noa.3  sample.1  sample.2  sample.3           noa        sample     ratio
          m*  445141.0  256474.0  142811.0  254850.0  107060.0  258338.0  89401.143052  75644.325430  0.846123
          mS   91268.0   62794.0   35692.0   54993.0   20419.0   50959.0  20378.661544  15329.317262  0.752224
           m   72454.0   49090.0   33596.0   52640.0   17860.0   34748.0  16741.602393  12675.370806  0.757118
         ..        ...       ...       ...       ...       ...       ...           ...           ...       ...
         WMM       NaN       2.0       2.0       1.0       NaN       8.0      0.741297      1.658007  2.236630
         WMW       NaN       1.0       2.0       2.0       NaN       2.0      0.557864      0.698816  1.252663
         MWW       NaN       NaN       2.0       NaN       4.0       2.0      0.748862      1.261906  1.685098
         [8842 rows x 9 columns]
        """
        counts = self.get_counts(
            pos=pos, how=how, **kwargs
        )  # compute position counts
        if exclude_empty:
            counts = counts[~counts.index.str.fullmatch(' *')]
        norm = counts.div(counts.sum()).mul(factor)  # normalize
        avg = norm.mean(axis=1)
        out = counts.assign(**{self.name: avg})
        if self.reference:
            ref = (
                self.reference.get_counts_ratio(
                    pos=pos, how=how, factor=factor, **kwargs
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
    def get_counts_ratio_pos(self, pos=None, how='aax', **kwargs):
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

        Examples
        --------

        Calculate the enrichement relative to the reference for the default -2/E/P/A positions.
         >>> sample.get_counts_ratio_pos()
         amino-acid         H         R         K  ...         W         *         m
         site                                      ...
         -2          1.062831  1.066174  1.012982  ...  1.046303       NaN  0.907140
         E           1.037079  1.018643  0.941939  ...  1.041217       NaN  0.933880
         P           1.093492  1.100380  1.045145  ...  1.107238       NaN  0.793043
         A           0.831129  1.005783  0.967491  ...  0.995833  1.143702  0.757118
         [4 rows x 22 columns]

        Calculate the enrichement relative to the reference for custom positions.
         >>> sample.get_counts_ratio_pos(('-3', '-2', 'E', 'P', 'A'))
         amino-acid         H         R         K  ...         W         *         m
         site                                      ...
         -3          1.032528  1.014771  0.987577  ...  1.045751       NaN  0.903142
         -2          1.062861  1.064912  1.013528  ...  1.048815       NaN  0.907531
         E           1.036543  1.018309  0.940488  ...  1.043321       NaN  0.934014
         P           1.092992  1.101174  1.045651  ...  1.106943       NaN  0.792341
         A           0.830804  1.005100  0.968136  ...  0.993697  1.143881  0.753449
         [5 rows x 22 columns]
        """

        if not pos:
            if how == 'nuc':
                pos = tuple(range(-6, 6))   # -2/E/P/A
            else:
                pos = ('-2', 'E', 'P', 'A')

        if isinstance(pos, str):
            pos = ranges(pos)

        aa_mode = how.startswith('aa')

        return (
            pd.concat(
                {
                    p: self.get_counts_ratio(p, how=how, **kwargs)['ratio']
                    for p in pos
                },
                axis=1,
            )
            .pipe(lambda x: x.reindex(aa_order) if aa_mode else x)
            .rename_axis(
                index='amino-acid' if aa_mode else 'nucleotide', columns='site'
            )
            .T
        )

    @lru_cache
    def DE(
        self,
        pos=None,
        how='aax',
        join=False,
        quiet=True,
        filter_size=True,
        multi=True,
        n_cpus=None,
        raw=False,
        _nocache=False,
        **kwargs,
    ):
        """
        Computes the differential expression between the sample and its reference.

        Parameters
        ----------
        pos : str
            Ribosome positions to consider to compute the differential expression.
        how : str, optional
            Type of inverse toeprints to analyze (see :meth:`Replicate.load_data`).
        join : bool, optional
            If True, joins the DE results back to the original `df`. Defaults to False.
        quiet : bool, optional
            If True, suppresses the console output of the `pydeseq2` library. Defaults to True.
        multi : bool, optional
            Whether to compute DE with a specific contrast (`cond` vs. `ref`). Defaults to True.
        n_cpus : int, optional
            The number of CPUs to utilize for parallel processing. Defaults to the total number of available CPUs.
        filter_size: bool
            Only considers reads for which an amino acid is present in all target positions.
        **kwargs: optional
            Additional parameters to :meth:`get_counts_ratio` and :meth:`Replicate.load_data`.
            For instance ``min_peptide`` and ``max_peptide`` are useful to filter the peptide size of the inverse toeprints to consider.

        Returns
        -------
        DataFrame
            DataFrame of the differential expression statistics with a row per motif.

        See Also
        --------
        Sample.get_counts_ratio: Gets the inverse toeprint counts and sample/reference ration of normalized counts.
        Sample.volcano: Draws a volcano plot from the Differential Expression data.
        Sample.subset_logo: Creates a logo from a subset of the Differential Expression data.

        Examples
        --------

        Compute the differential expression for positions E-P
         >>> sample.DE('E:P')
                 baseMean  log2FoldChange     lfcSE       stat        pvalue          padj  log10pvalue  log10padj
         QK   5537.704183        0.778031  0.073280  10.617238  2.477833e-26  7.582170e-24    25.605928  23.120206
         VI   6874.891363        0.295160  0.371018   0.795542  4.262985e-01  7.718778e-01     0.370286   0.112451
         MY    747.538317        0.263705  0.074294   3.549477  3.859965e-04           NaN     3.413417        NaN
         YY   2216.501684        0.259860  0.068213   3.809545  1.392226e-04  6.086018e-03     3.856290   2.215667
         WM    200.446070        0.226720  0.111555   2.032371  4.211614e-02           NaN     1.375551        NaN
         ..           ...             ...       ...        ...           ...           ...          ...        ...
         TP  15256.234795       -0.255940  0.061793  -4.141886  3.444618e-05  2.635133e-03     4.462859   2.579197
         mK  10824.395771       -0.308538  0.210353  -1.466765  1.424400e-01  4.737680e-01     0.846368   0.324434
         EP   8950.363473       -0.321266  0.068514  -4.689045  2.744828e-06  2.799725e-04     5.561485   3.552885
         PP  20880.530910       -0.372203  0.078851  -4.720365  2.354220e-06  2.799725e-04     5.628153   3.552885
         KK   7645.111411       -0.390365  0.096140  -4.060381  4.899280e-05  2.998359e-03     4.309868   2.523116
         [420 rows x 8 columns]
        """
        print(f'DE: {self.name}, {pos=}')

        if not self.reference:
            raise ValueError('cannot compute DE without a reference')

        if len(self.replicates) < 2:
            print(f'DE: {self.name} has not enough replicates to compute DE.')
            return None
            # pd.DataFrame(columns=['log2FoldChange', 'log10pvalue', 'log10padj'])

        cond = self.name
        ref = self.reference.name

        # df = self.reference.get_counts(pos=pos, how=how, **kwargs).join(
        #    self.get_counts(pos=pos, how=how, **kwargs)
        # )
        df = self.get_counts_ratio(pos, how=how, **kwargs)

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

        res = DE(
            df,
            sample_df,
            cond=cond,
            ref=ref,
            quiet=quiet,
            multi=multi,
            n_cpus=n_cpus,
            raw=raw,
            _nocache=_nocache,
            _cache_dir=self._cache_dir,
            _cache_prefix=_cache_prefix,
        )
        if join:
            return df.join(res, how='right')
        return res

    def hmap(
        self,
        r=None,
        c=None,
        *,
        pos=None,
        # how='aax',
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

        Examples
        --------
        Create a heatmap for positions E-P-A:

         >>> sample.hmap('E', 'A')

        .. image:: /_static/sample_hmap.png
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
        table = table.loc[~table.index.str.contains(' ')]

        names = list(map(get_ribosome_site, names))

        table.index = pd.MultiIndex.from_tuples(
            [tuple(x) for x in table.index], names=names
        )

        if not ax:
            _, ax = plt.subplots()

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
        # how='aax',
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

        Examples
        --------
        Create the default heatmap grid for all combinations of -2/E/P/A:

         >>> sample.hmap_grid()

        .. image:: /_static/sample_hmap_grid.png

        Create a heatmap grid for combinations of E/P/A:

         >>> sample.hmap_grid(['E', 'P', 'A'])

        .. image:: /_static/sample_hmap_grid_EPA.png
        """

        if not pos:
            pos = ['-2', 'E', 'P', 'A']

        grid_size = len(pos) - 1
        f, axes = plt.subplots(
            nrows=grid_size,
            ncols=grid_size,
            figsize=(grid_size * 5, grid_size * 5),
        )

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
                    # how=how,
                    col=col,
                    transform=transform,
                    cmap=cmap,
                    vmax=vmax,
                    center=center,
                    ax=ax,
                    **kwargs,
                )
        return f

    def hmap_pos(
        self,
        pos=None,
        *,
        # how='aax',
        # col='auto',
        # transform=np.log2,
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
        """
        Draws a volcano plot from the Differential Expression data.

        Parameters
        ----------
        pos : str, optional
            Positions used to compute the :meth:`DE`.
        query : str, optional
            Query used to select points to highlight and annotate by data.
        motif : str, optional
            Query used to select points to highlight by motif.
        ax : matplotlib.axes.Axes, optional
            ax to use for plotting, otherwise create a new figure.
        x : str, optional
            Column from the :meth:`DE` output to use as the x-axis.
        y : str, optional
            Column from the :meth:`DE` output to use as the y-axis.
        query_color : str, optional
            Color of the points in the query.
        motif_color : str, optional
            Color of the points in the motif.
        color : str, optional
            Color of the background points.
        params : dict, optional
            Optional parameters passed to scatter.
        annotate : bool, optional
            If True, annotate the points from the query.
        text_stroke : bool, optional
        outfile : str, optional
            If specified, save the figure to a file.
        density_thresh : int, optional
            Filters out annotation in dense parts of the plot. A low value decreases the number of annotations.
        kwargs : optional
            Optional parameters passed to :meth:`DE`.

        Returns
        -------
        matplotlib.axes.Axes

        See Also
        --------
        Sample.DE: Computes the differential expression between the sample and its reference.
        Sample.subset_logo: Creates a logo from a subset of the Differential Expression data.

        Examples
        --------
        Create a volcano plot for positions E-A:

         >>> sample.volcano('E:A')

        .. image:: /_static/sample_volcano.png

        Create a volcano plot for positions -2/E/P/A.
        Highlight and annotate the points above a threshold abs(L2FC) and log10(p-value):

         >>> sample.volcano('-2:A', query='(abs(log2FoldChange) > 2) & (log10pvalue > 10)')

        .. image:: /_static/sample_volcano_query.png

        Create a volcano plot for positions -2/E/P/A.
        Highlight the points for motifs containing a central QK motif.

         >>> sample.volcano('-2:A', motif='.QK.', motif_color='#FFAC1E',
                            query='', annotate=False)

        .. image:: /_static/sample_volcano_motif.png
        """
        df = self.DE(pos=pos, **kwargs)
        if df is None:
            print(f'Skipping volcano plot for {self}')
            return None
        # df.astype(float)  ### FIXME

        df['rank'] = df[x].abs() * df[y]
        df['qrank'] = df['rank'].rank(pct=True)
        df.sort_values(by='qrank', ascending=False, inplace=True)

        if not ax:
            _, ax = plt.subplots()

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
            query = 'keep'

        elif query.startswith('auto'):
            quantile = 0.90
            m = re.search(r'auto:\s*(\d+\.?\d*)', query)
            if m:
                quantile = float(m.group(1))
                query = 'auto'
            if query == 'auto':
                x_thresh = df[x].abs().quantile(quantile)
                y_thresh = df[y].abs().quantile(quantile)

                x_min = df[x].agg(['min', 'max']).abs().min()
                x_thresh = min(x_thresh, x_min)

                query = f'(log10pvalue > {y_thresh}) & (abs(log2FoldChange) > {x_thresh})'
            else:
                thresh = {2: (0.05, 0.2), 3: (0.05, 1), 4: (0.2, 2)}
                padj, l2fc = thresh.get(len(df.index[0]), (0.05, 1))
                query = f'(padj < {padj}) & (log2FoldChange > {l2fc})'

        if query:
            df_query = df.query(
                query, engine='python'
            )  # engine='python' to prevent bug with <NA> type
            df_query.plot.scatter(x=x, y=y, color=query_color, ax=ax, **params)
            if (annotate is True) or ('query' in annotate):
                dfs.append(df_query)
                colors.append(query_color)

        if motif:
            df_motif = df[df.index.str.match(motif)]
            df_motif.plot.scatter(x=x, y=y, color=motif_color, ax=ax, **params)
            if (annotate is True) or ('motif' in annotate):
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
                    left = ((sx < 0) and (sy < -sx)) or (0 < sx < sy)

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

    def subset_logo(
        self,
        pos,
        *,
        how='aax',
        query='',
        motif=None,
        logo_type='extra_counts_bits',
        ax=None,
        logo_kwargs=None,
        return_matrix=False,
        **kwargs,
    ):
        """
        Creates a logo from a subset of the Differential Expression data.

        This first runs the :meth:`.DE` method with ``pos``, then uses the ``query`` and ``motif`` parameters to filter
        the output, and uses the filtered table to produce a logo.
        Each motif is weighted using the average normalized counts from :meth:`.get_counts_ratio`.

        Parameters
        ----------

        pos: str
            Ribosome positions to consider to compute the differential expression. Passed to :meth:`DE`.
        how: str, optional
            Type of inverse toeprints to analyze (see :meth:`Replicate.load_data`).
        query: str, optional
            Query used to select rows in the :meth:`DE` output.
        motif: regex, optional
            Regex used to select motifs (e.g., '.QK.' in a 4 amino-acid motif would fix the central QK).
        logo_type: str, optional
            Type of logo to compute:

            - "raw_freq": unweighted frequencies of the amino-acids for all present motifs
            - "extra_counts": Computes the sum of extra counts (sample - reference) for each residue per position
            - "sum_log2FC": sum of the log2FoldChange for each residue per position
            - "<logo_type>_bits": If any of the above has a "_bits" suffix, an extra conversion to bits is performed.

        ax: matplotlib Axes, optional
            If passed, the figure will be drawn on the given Axes. A new Axes is created otherwise.
        logo_kwargs: dict, optional
            Additional parameters passed to logomaker.Logo.
        return_matrix: bool, optional
            If True, the logo matrix is returned as together with the logo as (logo, matrix).
        kwargs: optional
            Additional parameters passed to the :meth:`.DE` method.

        Returns
        -------
        logomaker.Logo

        Examples
        --------
         >>> sample.subset_logo('-2:A', query='(log2FoldChange > 2) & (log10pvalue > 2)')

        .. image:: /_static/sample_subset_logo.png
        """

        if logo_kwargs is None:
            logo_kwargs = {
                'color_scheme': aa_colors
                if how.startswith('aa')
                else 'classic'
            }

        df = self.DE(pos, how=how, join=True, **kwargs)

        if motif is not None:
            df = df[df.index.str.match(motif)]

        if query:
            df = df.query(query)

        return motif_logo(
            df,
            ref_spl=[self.reference.name, self.name],
            ax=ax,
            logo_type=logo_type,
            logo_kwargs=logo_kwargs,
            return_matrix=return_matrix,
        )

    def all_logos(self, logo_kwargs=None, **kwargs):
        """
        Creates a logo for all positions for each replicate in the sample.

        Parameters
        ----------
        logo_kwargs : dict, optional
            Optional parameters to pass to logomaker.Logo.
        kwargs : optional
            Optional parameters to pass to :meth:`logo`.

        See Also
        --------
        Sample.logo: Creates a logo for the selected positions.

        Examples
        --------
         >>> tcx2.all_logos()

        .. image:: /_static/sample_all_logos.png
        """
        _, axes = plt.subplots(nrows=len(self.replicates), sharex=True)
        for r, ax in zip(self.replicates, axes.flat):
            r.logo(logo_kwargs=logo_kwargs, ax=ax, **kwargs)
            ax.set_title(r)

    def logo(
        self,
        pos=None,
        logo_kwargs=None,
        ax=None,
        vs_ref=True,
        how='aax',
        **kwargs,
    ):
        """
        Creates a logo for the selected positions.

        Parameters
        ----------
        pos : tuple
            Positions to use in the logo (see :meth:`get_counts_ratio_pos`)
        logo_kwargs : dict, optional
            Optional parameters to pass to logomaker.Logo.
        ax : matplotlib.axes.Axes, optional
            ax to use for plotting, otherwise create a new figure.
        kwargs : optional
            Optional parameters to pass to :meth:`get_counts_ratio_pos`.

        See Also
        --------
        Sample.all_logos: Creates a logo for all positions for each replicate in the sample.
        Sample.subset_logo: Creates a logo from a subset of the Differential Expression data.

        Examples
        --------
         >>> tcx2.logo()

        .. image:: /_static/sample_logo.png
        """
        if vs_ref:
            df = np.log2(
                self.get_counts_ratio_pos(pos=pos, how=how, **kwargs)
            ).fillna(0)
        else:
            raise NotImplementedError('vs_ref=False is not implemented yet')

        import logomaker

        if logo_kwargs is None:
            logo_kwargs = {
                'color_scheme': 'NajafabadiEtAl2017'
                if how.startswith('aa')
                else 'classic',
                'flip_below': False,
            }

        # df = logomaker.transform_matrix(df, from_type='counts', to_type=type)
        # print(df)

        if not ax:
            ax = plt.subplot()

        logo = logomaker.Logo(df.reset_index(drop=True), ax=ax, **logo_kwargs)
        logo.ax.set_xticks(range(df.shape[0]), df.index)
        return logo

    @cached_property   # FIXME consider not using a property if we want to be able to modify replicates
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

        Examples
        --------
         >>> sample.itp_len
              length replicate     count sample
         0        51         1  115732.0    spl
         1        20         1  444506.0    spl
         2        41         1  130495.0    spl
         3        23         1  198257.0    spl
         4        17         1   55786.0    spl
         ..      ...       ...       ...    ...
         328     106         3       NaN    spl
         329     143         3       NaN    spl
         330     102         3       NaN    spl
         331     104         3       NaN    spl
         332     221         3       NaN    spl
         [333 rows x 4 columns]
        """

        return (
            pd.concat(
                {
                    r.replicate: pd.Series(r.infos()['lengths'])
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

        See Also
        --------
        DataSet.itp_len_plot

        Notes
        -----
        - The x-axis represents the distance from the 3' end of the inverse-toeprint in nucleotides.
        - The y-axis shows the counts of inverse-toeprints, either absolute or normalized per million reads.
        - Each replicate is plotted independently and distinguished by the `hue` attribute in the plot.

        Examples
        --------
         >>> sample.itp_len_plot()

        .. image:: /_static/sample_itp_len_plot.png
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
            _, ax = plt.subplots()

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
        """DataFrame of the counts of each inverse toeprint length per Replicate"""
        return (
            pd.concat(
                {
                    r.name: pd.Series(r.infos()['lengths'])
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
        """
        Plots a virtual inverse-toeprint gel.

        Parameters
        ----------
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
        interactive : bool, optional
            If True, creates an interactive display to test all parameters.
            Requires to run in a notebook with ipywidgets.

        Returns
        -------
        matplotlib.axes.Axes

        Examples
        --------
        Create a virtual inverse-toeprint with bands.

         >>> sample.itoeprint(exposure=3)

        .. image:: /_static/sample_itoeprint.png

        Create a virtual inverse-toeprint with shades:

         >>> sample.itoeprint(plot='shades')

        .. image:: /_static/sample_itoeprint_shades.png

        Create an interactive visualization (requires running in a notebook with ipywidgets):

         >>> sample.itoeprint(interactive=True)

        .. image:: /_static/sample_itoeprint_interactive.png
        """
        if interactive:
            import ipywidgets as widgets
            from ipywidgets import fixed, interact

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
            return None
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
    __doc__ = (
        rf"""
    Loads an iTP-Seq dataset and provides methods for analyzing and visualizing the data.

    A DataSet object is constructed to handle iTP-Seq Samples with their respective Replicates.
    By default, it infers the files to uses in the provided directory by looking for "\*.{ITP_FILE_SUFFIX}.json" files produced
    during the initial step of pre-processing and filtering the fastq files.
    It uses the pattern of the file names to group the Replicates into a Sample, and to define which condition
    is the reference in the DataSet (the Sample with name "noa" by default).

    Attributes
    ----------
    data_path : str or Path
        Path to the data directory containing the output files from the fastq pre-processing.
    result_path: str or Path
        Path to the directory where the results of the analysis will be saved.
    samples: list or dict or None
        List or dictionary of Samples in the DataSet.
        By default, it is None and will be populated automatically if data_path is provided.
    keys: tuple
        Properties in the file name to use for identifying the reference.
    ref_labels: str or tuple
        Specifies the reference: e.g. 'noa' or (('sample', 'noa'),)
    cache_path: str or Path
        Path used to cache intermediate results. By default, this creates a subdirectory called "cache" in the result_path directory.
    file_pattern: str
        Regex pattern used to identify the sample files in the data_path directory.
        If None, defaults to `r'(?P<lib_type>[^_]+)_(?P<sample>[^_\d]+)(?P<replicate>\d+)'`
        which matches files like nnn15_noa1.{ITP_FILE_SUFFIX}.json, nnn15_tcx2.{ITP_FILE_SUFFIX}.json, etc.
    allow_partial_keys: bool
        If no exact match of the keys if found, try to map a reference using partial keys.
    ref_mapping: dict or None
        If set, do not try to infer the references but use the passed dictionary as mapping.
        The dictionary should have a format: `{{'sample.id': 'ref.id'}}` where "sample.id" and "ref.id" are the labels generated upon import.
    """
        + """
    Examples
    --------
    Creating a DataSet from a simple antibiotic treatment (tcx) vs no treatement (noa) with 3 replicates each (1, 2, 3).

    Load a dataset from the current directory, inferring the samples automatically.
     >>> from itpseq import DataSet
     >>> data = DataSet(data_path='.')
     >>> data
     DataSet(data_path=PosixPath('.'),
             file_pattern='(?P<lib_type>[^_]+)_(?P<sample>[^_\\d]+)(?P<replicate>\\d+)',
             samples=[Sample(nnn15.noa:[1, 2, 3]),
                      Sample(nnn15.tcx:[1, 2, 3], ref: nnn15.noa)],
             )

    Same as above, but only use "sample" as key.
     >>> data = DataSet(data_path='.', keys=['sample'])
     >>> data
     DataSet(data_path=PosixPath('.'),
             file_pattern='(?P<lib_type>[^_]+)_(?P<sample>[^_\\d]+)(?P<replicate>\\d+)',
             samples=[Sample(noa:[1, 2, 3]),
                      Sample(tcx:[1, 2, 3], ref: noa)],
             )

    Compute a standard report and export it as PDF
     >>> data.report('my_experiment.pdf')

    Display a graph of the inverse-toeprints lengths for each sample
     >>> data.itp_len_plot(row='sample')
    """
    )

    def __init__(
        self,
        data=None,
        *,
        data_path: Path = None,
        result_path: Path = None,
        samples: Optional[dict] = None,
        keys=None,
        ref_labels: Optional[Union[str, tuple]] = 'noa',
        # max_workers: int = 4,
        cache_path=None,
        file_pattern=None,
        allow_partial_keys=True,
        ref_mapping=None,
    ):
        # if data is passed, try to infer if this is data_path or samples:
        if data:
            if data_path or samples:
                raise ValueError(
                    'If data is set, cannot use data_path or samples'
                )
            if isinstance(data, (str, Path)):
                data_path = data
            elif isinstance(data, (tuple, list, dict)):
                samples = data
            else:
                raise ValueError(f'Unrecognized data: {type(data)}')

        self.data_path = Path(data_path) if data_path else None
        if result_path is None and self.data_path:
            result_path = self.data_path / 'results'
        self.result_path = Path(result_path) if result_path else None
        if cache_path is None:
            if self.result_path:
                cache_path = self.result_path / 'cache'
            else:
                import tempfile

                cache_path = tempfile.TemporaryDirectory().name
                print(f'Creating temporary cache directory: "{cache_path}"')
        if self.result_path and not self.result_path.exists():
            self.result_path.mkdir(parents=True)
        self.cache_path = Path(cache_path)
        if not self.cache_path.exists():
            self.cache_path.mkdir(parents=True)
        # self.file_pattern = r'nnn15_(?P<sample>[^_]+)(?P<replicate>\d+)'
        self.file_pattern = file_pattern or FILE_PATTERN
        # self.aafile_pattern = (
        #     aafile_pattern or f'{lib_type}_{sample}{replicate}.aa.{ITP_FILE_SUFFIX}.txt'
        # )
        # if no keys are provided, use all the named capturing groups of file_pattern, excepted "replicate"
        self.keys = (
            keys
            if keys is not None
            else [
                g
                for g in re.compile(self.file_pattern).groupindex
                if g != 'replicate'
            ]
        )
        self.ref_mapping = ref_mapping

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

        # Define the Samples in the DataSet
        # if data_path is provided, try to infer the samples from the file names using file_pattern
        if data_path:
            if not samples:
                self.samples = self._infer_samples(allow_partial_keys)
            else:
                raise ValueError('Cannot use both "data_path" and "samples".')
        # if a list or dictionary of Samples/Replicates is provided
        # create the objects if needed and assign them to the DataSet
        elif samples:
            if isinstance(samples, list):
                samples = {
                    f'sample{i}': S for i, S in enumerate(samples, start=1)
                }

            self.samples = {
                k: data
                if isinstance(data, Sample)
                else Sample(replicates=data, name=k, dataset=self)
                for k, data in samples.items()
            }
        else:
            self.samples = {}

        if self.ref_mapping is not None:
            self.set_references(self.ref_mapping)

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
        data_path = (
            f"data_path='{self.data_path}'{filepat}" if self.data_path else ''
        )
        samples = (
            f'{sep if self.data_path else ""}samples=[{sep2.join(repr(s) for s in self.samples.values())}]{sep}'
            if self.samples
            else ''
        )
        # return f'''DataSet(data_path={repr(self.data_path)}, result_path={repr(self.result_path)}{samples})'''
        return f"""DataSet({data_path}{samples})"""

    def _clear_cache(self, force=False):
        import os
        from shutil import rmtree

        if self.cache_path.exists() and (
            force
            or input(
                f'Delete {len(next(os.walk(self.cache_path))[2])} files in "{self.cache_path}"? (y/N): '
            ).lower()
            == 'y'
        ):
            rmtree(self.cache_path)

    @property
    def samples_with_ref(self):
        """Dictionary of the samples that have a reference"""
        return {k: s for k, s in self.samples.items() if s.reference}

    def set_references(self, ref_mapping=None, exact_mapping=False):
        """
        Sets the Sample references from a mapping

        Parameters
        ----------
        ref_mapping : mapping, optional
            Mapping (ex. dictionary) of {'sample_label': Replicate or 'replicate_label'}
        exact_mapping : bool, optional
            If set to True, will remove the existing references for samples not defined in the mapping
        """
        if ref_mapping is not None:
            for k, s in self.samples.items():
                if k in ref_mapping:
                    ref = ref_mapping[k]
                    if not isinstance(ref, Replicate):
                        ref = self.samples[ref]
                    s.reference = ref
                elif exact_mapping:
                    s.reference = None

    def _infer_samples(self, allow_partial_keys=True):
        """Infers sample names from the files in the data path."""
        inferred_samples = defaultdict(list)
        # file_paths = list(self.data_path.glob("*.json"))

        for f in sorted(
            self.data_path.iterdir()
        ):  ## TODO: wrap with SampleGrouper
            if m := re.search(
                self.file_pattern + rf'(?=\.{ITP_FILE_SUFFIX}\.json$)', f.name
            ):
                labels = m.groupdict()
                # if set(labels) >= {'sample', 'replicate'}:
                if set(labels) > {'replicate'}:
                    rep = labels['replicate']
                    # keys = tuple(sorted(t for t in labels.items() if t[0] in self.keys))
                    keys = dict_to_tuple(labels, keep=self.keys)
                    inferred_samples[keys].append(
                        {
                            #'file_prefix': str(self.data_path / self.aafile_pattern.format(**labels)),
                            'file_prefix': str(self.data_path / m.group(0)),
                            'replicate': rep,
                            'labels': labels,
                        }
                    )

        # initialize Samples (without reference)
        samples = {}
        for key, data in inferred_samples.items():
            s = Sample(
                labels=dict(key), replicates=data, dataset=self, keys=self.keys
            )
            samples[s.name] = s

        if self.ref_mapping is None:
            # create a dictionary of reference samples based on self.ref_labels
            ref_keys = set(dict(self.ref_labels))
            ref_items = dict(self.ref_labels).items()
            ref_samples = {
                dict_to_tuple(s.labels, ignore=ref_keys): s
                for k, s in list(samples.items())[::-1]
                if ref_items <= s.labels.items()
            }

            # create a dictionary of reference samples
            # ignoring keys with None as value
            # this is used if no exact match is found
            ref_samples_minkey = defaultdict(list)
            for k, s in ref_samples.items():
                new_key = tuple(t for t in k if t[0] in ref_keys)
                ref_samples_minkey[new_key].append(s)

            # ensure there is a single value for a given key
            # if several references match a single key, they will be removed from ref_samples_minkey
            for k in list(ref_samples_minkey):
                if len(ref_samples_minkey[k]) > 1:
                    # spl_str = ', '.join(str(x) for x in ref_samples_minkey[k])
                    # print(f'Multiple references for {dict(k)}: [{spl_str}]')
                    ref_samples_minkey.pop(k)
            # flatten dictionary to keep the single references
            ref_samples_minkey = {
                k: l[0] for k, l in ref_samples_minkey.items()
            }

            # loop over the samples
            # first try to find a reference with an exact match of the keys
            # if no exact match is found, try to find a reference with less keys that all match
            ref_ids = {s.name for s in ref_samples.values()}
            for k, s in samples.items():
                if k in ref_ids:  # if this is a reference sample, continue
                    continue
                # items = dict(dict_to_tuple(s.labels, ignore=ref_keys)).items()
                sample_labels = dict_to_tuple(s.labels, ignore=ref_keys)

                # exact match
                if sample_labels in ref_samples:
                    # print(f'{s} -> {ref_samples[sample_labels]}')
                    s.reference = ref_samples[sample_labels]
                    continue
                spl_lbl = dict(sample_labels).items()

                # sub match (the reference can have less keys than the sample, but all must match)
                if allow_partial_keys:
                    for ref_lbl, ref_s in ref_samples_minkey.items():
                        ref_lbl = dict(ref_lbl).items()
                        if ref_lbl <= spl_lbl:
                            # print(f'{s} -> {ref_s}')
                            s.reference = ref_s
                            continue

        return {k: samples[k] for k in sorted(samples)}

    def reorder_samples(self, order, validate=True, reorder_replicates=True):
        """
        Reorders the samples in the DataSet.

        This method is useful to specify a custom order to use in the different graphs
        (e.g. in :meth:`itoeprint`).

        Parameters
        ----------
        order : list
            New order of the samples
        validate : bool
            Ensure the new order contains all samples.
        reorder_replicates : bool
            Also reorder the replicates.

        Examples
        --------
         >>> data.reorder_samples(['sampleY', 'sampleX', 'reference'])
        """
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
        """
        Displays summary information about the dataset NGS reads per replicate.

        This information is computed during the parsing step and includes:

        - the total number of reads,
        - the number of reads without adaptors,
        - the number of reads that are contaminants,
        - the number of reads with a low quality,
        - the number of reads that are too short or too long,
        - the number of extra nucleotides at the 3'-end of the inverse-toeprints,

        Parameters
        ----------
        html : bool
            if True, returns the table as HTML, otherwise as DataFrame (default).

        Examples
        --------
         >>> dataset.infos()
                   total_sequences  noadaptor  contaminant  lowqual  tooshort  toolong   extra0   extra1   extra2  MAX_LEN
         noa.1             9036255     799955         1219   700374    299502  3376581  2434092  2709762  3092446       44
         noa.2             8154560     407750         1318   680813    154587  4158921  2329190  2582052  2835568       44
         noa.3             7725561     623037         1065   353401    279104  3505909  2216460  2402957  2483107       44
         sample.1          8384889     714414         1192   685017    385537  3341987  2308291  2528638  2833546       44
         sample.2          9120203     498202         1659   513308    104071  5664107  2673062  2972850  2976089       44
         sample.3          8490958    1043590         1328   409697    187746  4004073  2243720  2555783  2647865       44
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
                {name: r.infos() for name, r in self.replicates.items()},
                orient='index',
            )
            .sort_index(axis=1, key=lambda x: x.map(key))
            .select_dtypes('number')
        )
        if html:
            return table_to_html(out)
        return out

    def DE(self, pos='E:A', **kwargs):
        """
        Computes the log2-FoldChange for each motif described by `pos` for each sample in the DataSet relative to their reference

        Attributes
        ----------
        pos: str
            position of the motif to consider.
            This ca be a range of positions (e.g. '-2:A' for -2/E/P/A sites)
            or a combination of disjoint positions (e.g. 'E,A' for the combination of E and A sites).
        kwargs:
            parameters passed to `Sample.get_counts` computes the counts for the given motif.
            for example `min_peptide=3` to consider only peptides of at least 3 amino acids.

        Returns
        -------
        dictionary of {sample_name: DE}

        See Also
        --------
        Sample.DE: Compute differential Expression for a Sample
        """
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
        limit=100,
        norm=True,
        hue='auto',
        plt_kwargs=dict(kind='line', height=2, aspect=3),
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
        col: string, optional
            attribute to use as columns in the FacetGrid
        row: string, optional
            attribute to use as rows in the FacetGrid
        min_codon : int, optional
            The minimum codon position to annotate on the plot. Defaults to 0.
        max_codon : int, optional
            The maximum codon position to annotate on the plot. Defaults to 10.
        limit : int, optional
            The maximum length to include in the plot. Defaults to 100.
        norm : bool, optional
            Whether to normalize counts to reads per million. Defaults to False.
        hue: str, optional
            Parameter to use a hue in the FacetGrid (by default 'replicate').
        plt_kwargs: dict, optional
            parameters used in the FacetGrid if col/row is used.

        Returns
        -------
        matplotlib.axes.Axes or seaborn.axisgrid.FacetGrid

        See Also
        --------
        Sample.itp_len_plot

        Notes
        -----
        - The x-axis represents the distance from the 3' end of the inverse-toeprint in nucleotides.
        - The y-axis shows the counts of inverse-toeprints, either absolute or normalized per million reads.
        - Each replicate is plotted independently and distinguished by the `hue` attribute in the plot.

        Examples
        --------
        Plot a line with error band for each sample:

         >>> dataset.itp_len_plot()

        .. image:: /_static/dataset_itp_len_plot.png

        Create a figure with a subplot per sample and a line per replicate:

         >>> dataset.itp_len_plot(row='sample')

        .. image:: /_static/dataset_itp_len_plot_row.png
        """

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
            for ax_ in g.axes.flat:
                itp_len_add_positions(
                    ax_, min_codon=min_codon, max_codon=max_codon
                )
            return g

        if not ax:
            _, ax = plt.subplots()

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
        """DataFrame of the counts of each inverse toeprint length per Replicate"""
        return (
            pd.concat(
                {
                    r.name: pd.Series(r.infos()['lengths'])
                    for _, r in self.replicates.items()
                },
                axis=1,
            )
            .rename(int)
            .sort_index(ascending=False)
        )

    itoeprint = Sample.itoeprint

    def report(self, template='report', output=None):
        """
        Create a report for the DataSet.

        Parameters
        ----------
        template : str, optional
            Which template to use for the report.
            There is currently only one template, but you can add custom ones in `itpseq/templates/`.
        output : str, optional
            Name of the output file in which to write the report.
            The type of the file in inferred from the extention.
            Only PDF (`.pdf)` and HTML (`.html`) are supported.

        Examples
        --------
        Create a PDF report
         >>> dataset.report(output='my_report.pdf')

        Create an HTML report
         >>> dataset.report(output='my_report.html')

        Create an HTML report in the notebook
         >>> dataset.report()
        """

        from jinja2 import Environment, FileSystemLoader

        template_dir = Path(__file__).resolve().parent / 'templates'
        static_dir = Path(__file__).resolve().parent / 'static'
        env = Environment(loader=FileSystemLoader(template_dir))
        if not template_dir.joinpath(f'{template}.html').exists():
            raise ValueError(
                f'Template {template}.html not found in {template_dir}'
            )
        html_template = env.get_template(f'{template}.html')

        # additional variables
        render_kwargs = {
            'today': datetime.date.today(),
        }

        if output:
            output = Path(output)
            if output.suffix == '.html':
                css_file = static_dir / f'{template}-html.css'
                if not css_file.exists():
                    css_file = static_dir / 'report-html.css'
                html = html_template.render(
                    dataset=self,
                    plot_to_html=plot_to_html,
                    css_file=css_file,
                    output='html',
                    **render_kwargs,
                )
                with open(output, 'w', encoding='utf-8') as f:
                    f.write(html)
            elif output.suffix == '.pdf':
                from weasyprint import CSS, HTML

                css_file = static_dir / f'{template}-pdf.css'
                if not css_file.exists():
                    css_file = static_dir / 'report-pdf.css'
                html = html_template.render(
                    dataset=self,
                    plot_to_html=plot_to_html,
                    output='pdf',
                    **render_kwargs,
                )
                HTML(string=html).write_pdf(
                    output, stylesheets=[CSS(filename=css_file)]
                )
            return None
        else:
            html = html_template.render(
                dataset=self,
                plot_to_html=plot_to_html,
                output='other',
                **render_kwargs,
            )

            if IPYTHON:
                from IPython.display import HTML

                return HTML(html)

            return html


if __name__ == '__main__':
    pass
