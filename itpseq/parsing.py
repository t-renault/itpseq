#!/usr/bin/env python3
"""Module to parse iTP-Seq fastq files for downstream analysis"""
import bz2
import gzip
import lzma
import sys
import zipfile
from contextlib import nullcontext
from functools import partial, wraps
from pathlib import Path

from .config import *

__all__ = ['parse', 'format_sequences']

DEFAULTS = dict(
    a1='GTATAAGGAGGAAAAAAT',
    a2='GGTATCTCGGTGTGACTG',
)

# fmt: off
codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}
# fmt: on

CONTAMINANTS = 'TCCAACATGCTGAGC|^GATCCTTTTTA'

#####

import numpy as np

to_aa = codon_table.copy()
to_aa['   '] = ' '
to_aa['\n'] = '\n'


def fastq_iterator(filename):
    """
    Reads 'filename' as a fastq sequence and yield only the sequence and quality lines

    This does not perform any check on the validity of the file!
    """
    from itertools import islice

    fn = str(filename)
    if fn.endswith('.gz'):
        import gzip

        opener = gzip.open
    elif fn.endswith('.bz2'):
        import bz2

        opener = bz2.open
    elif fn.endswith('.zip'):
        import zipfile

        opener = zipfile.ZipFile
    elif fn.endswith('.xz') or fn.endswith('.lzma'):
        import lzma

        opener = lzma.open
    else:
        opener = open

    with opener(filename, 'rt') as f:
        while True:
            chunk = list(islice(f, 4))
            if not chunk:
                break
            _, seq, _, score = chunk
            yield (seq[:-1], score[:-1])


def parse_filter_fastq(
    filename,
    *,
    a1=DEFAULTS['a1'],  # adaptators
    a2=DEFAULTS['a2'],
    mm1=2,  # allowed number of mismatches
    mm2=2,
    limit=None,  # max sequences to read
    min_seq_len=3 * 1,
    max_seq_len=3 * 10,
    quality=30,
    start='ATG',
    contaminants=CONTAMINANTS,
    untranslated_overhang=12,  # minimum extra nt after A-site (for RNAse R w/ ribosome this is 12-14)
    **kwargs,  # catch all other parameters
):
    """
    Takes a 'filename' as input (fastq format), loops over the fastq using
    fastq_iterator and performs several checks to extract valid ITP sequences.
    Also computes various statistics on the dataset.

    Parameters
    ----------

    filename:
        fastq file to process
    a1: str
        Sequence of the left adaptator
    a2: str
        Sequence of right adaptator
    mm1: int
        Number of tolerated mismatches in a1
    mm2: int
        Number of tolerated mismatches in a2
    limit: int or None
        Maximum number of sequences to process (useful for quick tests)
    min_seq_len: int or None
        Minimum length of the matched sequence to keep
    max_seq_len: int or None
        Maximum length of the matched sequence to keep
    quality: int
        Minimum quality required in the coding sequence to keep a read
    contaminants: str,
        Regex defining the contaminants to remove the matching reads
    untranslated_overhang: int
        Minimum extra nucleotides after the A-site (for RNAse R w/ ribosome this is 12).
        The possible overhangs are x / x+2.
    """

    from collections import Counter

    import regex

    items = [
        'total_sequences',
        'noadaptor',
        'lowqual',
        'contaminant',
        'tooshort',
        'toolong',
        'MAX_LEN',
    ]
    stats = Counter(dict.fromkeys(items, 0))   # initialize order of the items
    minquals = Counter()
    lengths = Counter()

    if not min_seq_len:
        min_seq_len = 0

    min_seq_len_overhang = min_seq_len + untranslated_overhang + 2

    if max_seq_len:
        max_seq_len_overhang = max_seq_len + untranslated_overhang + 2
    else:
        max_seq_len_overhang = float('inf')

    seqs = []
    MAX_LEN = 0
    extras = []

    def parse(record):
        keep = True
        seq, qual = record
        m = REG.search(seq)

        stats['total_sequences'] += 1

        if m:
            subseq = m.group()
            subqual = qual[slice(*m.span())]
            if max_seq_len:   # limit quality filter to first x codons
                subqual = subqual[:max_seq_len]
            minqual = ord(min(subqual)) - 33

            # span, subsequence, quality substring, min quality
            # out.append((*m.span(), subseq, subqual, minqual))
            subseq_len = len(subseq)
            lengths[subseq_len] += 1

            extra = subseq_len % 3
            MAX = subseq_len + (2 - extra)
            subseq = subseq + ' ' * (2 - extra)
            stats[f'extra{extra}'] += 1

            # As in first 18 nucleotides
            # if subseq[:18].count('A') >= 18:
            if subseq[:18] == 'AAAAAAAAAAAAAAAAAA':
                stats['A_stretch'] += 1
                keep = False

            if REG_conta.search(subseq):
                stats['contaminant'] += 1
                keep = False

            minquals[minqual] += 1

            if minqual < quality:
                stats['lowqual'] += 1
                keep = False
            if 'N' in subseq:
                stats['subseq_N'] += 1
                keep = False

            if subseq_len < min_seq_len_overhang:
                stats['tooshort'] += 1
                keep = False
            elif subseq_len > max_seq_len_overhang:
                stats['toolong'] += 1
                keep = False

            if keep:
                nonlocal MAX_LEN
                MAX_LEN = max(MAX, MAX_LEN)
                nonlocal extras
                extras.append(extra)
                return subseq

        else:
            stats['noadaptor'] += 1
        if 'N' in seq:
            stats['seq_N'] += 1

    records = fastq_iterator(filename)

    if limit:
        from itertools import islice

        records = islice(records, int(limit))

    REG = regex.compile(
        f'(?<=(?:{a1}){{s<={mm1}}}){start}.+(?=(?:{a2}){{s<={mm2}}})',
        regex.BESTMATCH,
    )
    REG_conta = regex.compile(f'(?:{contaminants}){{e<={2}}}')

    for record in records:
        seq = parse(record)
        if seq:
            seqs.append(seq)

    stats['MAX_LEN'] = MAX_LEN
    stats['min_quality_counts'] = minquals
    stats['lengths'] = lengths

    return seqs, stats, extras


def parse_file_wrapper(filename, **global_kwargs):
    """Wrapper to :func:`parse_filter_fastq` and :func:`export_data` for multiprocessing"""
    seqs, stats, extras = parse_filter_fastq(filename, **global_kwargs)
    if global_kwargs.get('save'):
        export_data(
            filename,
            seqs=seqs,
            stats=stats,
            outdir=global_kwargs.get('outdir'),
        )
    return seqs, stats, extras


def parse_all(files=None, pattern=None, save=False, outdir=None, **kwargs):
    """
    Apply parse_filter_fastq on multiple files
    """

    from multiprocessing import Pool, cpu_count

    if pattern:
        from glob import glob

        files = glob(pattern)

    f = partial(parse_file_wrapper, **kwargs, save=save, outdir=outdir)

    with Pool(max(1, cpu_count() - 2)) as pool:
        result = pool.map(f, files)

    short = [Path(f).name for f in files]

    return dict(zip(short, result))


def simple_graph(
    counter,
    vmin=None,
    vmax=None,
    start='auto',
    trim_below=0.02,
    ticks=10,
    wrap=80,
):
    """Generate an ASCII bar graph from a Counter"""
    import math

    bars = ' _▁▂▃▄▅▆▇█'[1:]
    max_x = max(counter)
    y_values = np.array([counter[i] for i in range(max_x)])
    max_y = y_values.max()

    # floating vmin/vmax = relative to max_y
    if isinstance(vmin, float):
        vmin = int(vmin * max_y)
    if isinstance(vmax, float):
        vmax = int(vmax * max_y)

    # trim below percent max_y
    if isinstance(trim_below, float):
        trim_below = int(trim_below * max_y)

    if trim_below:
        l = len(y_values)
        y_values[y_values < trim_below] = 0
        y_values = np.trim_zeros(y_values, trim='b')
        max_x -= l - len(y_values)

    if start == 'auto':
        start = np.min(np.nonzero(np.hstack((y_values, 1)))) // ticks * ticks
    max_x -= start
    y_values = y_values[start:]

    if not vmin:
        vmin = 0
    if not vmax:
        vmax = max_y
    # handle custom vmax
    if vmax < max_y:
        bars += '▒'
        vmax += 1

    b = (
        ((y_values - vmin) / (vmax - vmin) * (len(bars) - 1))
        .astype(int)
        .clip(0, len(bars) - 1)
        .tolist()
    )

    graph = ''.join(bars[i] for i in b)
    xticks = (f'╹{" "*(ticks-1)}' * math.ceil((max_x + start) / ticks))[:max_x]

    # left-aligned tick labels
    # labels = ''.join(f'{i: <{ticks}}' for i in range(start, L+start, ticks))

    # right-aligned tick labels
    labels = ''.join(
        f'{i: >{ticks}}' for i in range(start + ticks, max_x + start, ticks)
    )
    S = str(start)
    labels = S + labels[len(S) - 1 :]

    out = [graph, xticks, labels]

    if wrap:
        import textwrap
        from itertools import chain

        wrapper = textwrap.TextWrapper(width=wrap)
        out = list(chain.from_iterable(zip(*(wrapper.wrap(s) for s in out))))

    return '\n'.join(out)


def seq2aa(seq):
    """
    Translates a nucleotide sequence into an amino acid sequence

    We are not using BioPython to improve the speed of the module
    """
    return ''.join(
        [to_aa.get(seq[i : i + 3], '?') for i in range(0, len(seq), 3)]
    )


def export_data(
    filename,
    *,
    seqs=None,
    stats=None,
    outdir=None,
    MAX=None,
    untranslated_overhang=12,
):
    __doc__ = rf"""
    Exports the itpseq output files from the parsed data
    
    This creates the following files:
    - ``<file_prefix>.nuc.{ITP_FILE_SUFFIX}.txt``: inverse toeprints as nucleotides
    - ``<file_prefix>.aa.{ITP_FILE_SUFFIX}.txt``: inverse toeprints as amino acids
    - ``<file_prefix>.{ITP_FILE_SUFFIX}.json``: metadata as JSON
    - ``<file_prefix>.{ITP_FILE_SUFFIX}.log``: log file
    """

    print(f'exporting data from: {filename}')

    fname = Path(filename)
    if not fname.exists():
        raise ValueError(f'"{fname}" does not exist')
    if not fname.is_file():
        raise ValueError(f'"{fname}" is not a file')

    if outdir:
        outdir = Path(outdir)
        if not outdir.exists():
            outdir.mkdir()
    else:
        outdir = fname.parent

    import re

    base_fname = outdir / re.sub(
        r'(\.assembled)?\.f(ast)?q(\.(gz|bz2|zip|xz|lzma))?$', '', fname.name
    )
    f_log = base_fname.with_suffix(f'.{ITP_FILE_SUFFIX}.log')
    f_json = base_fname.with_suffix(f'.{ITP_FILE_SUFFIX}.json')
    # f_extra = base_fname.with_suffix(f'.{ITP_FILE_SUFFIX}.extra')
    f_seq_nuc = base_fname.with_suffix(f'.nuc.{ITP_FILE_SUFFIX}.txt')
    f_seq_aa = base_fname.with_suffix(f'.aa.{ITP_FILE_SUFFIX}.txt')

    max_untranslated_overhang = untranslated_overhang + 2

    if stats:
        with open(f_json, 'w', encoding='utf-8') as f_json:
            import json

            f_json.write(json.dumps(stats))

        with open(f_log, 'w', encoding='utf-8') as f_log:
            f_log.write('')
            N = len(str(stats['total_sequences'])) + 2
            f_log.write(
                f"Total sequences read: {stats['total_sequences']: >{N}}\n\n"
                f'Rejected sequences\n'
                f"    no adaptor found: {stats['noadaptor']: >{N}}\n"
                f"         contaminant: {stats['contaminant']: >{N}}\n"
                f"      stretches of A: {stats['A_stretch']: >{N}}\n"
                f"      seq contains N: {stats['seq_N']: >{N}}\n"
                f"   subseq contains N: {stats['subseq_N']: >{N}}\n"
                f"         low quality: {stats['lowqual']: >{N}}\n"
                f"    subseq too short: {stats['tooshort']: >{N}}\n"
                f"     subseq too long: {stats['toolong']: >{N}}\n\n"
            )
            if stats['lengths'].values():
                f_log.write(
                    'Distribution of subsequences length:\n'
                    f"scale: {min(stats['lengths'].values())}–{max(stats['lengths'].values())}\n"
                )
            else:
                f_log.write(
                    'Distribution of subsequences length:\nno sequences\n'
                )
            f_log.write(simple_graph(stats['lengths'], start=1) + '\n')
    # if extras:
    #    np.array(extras, dtype=np.int8).tofile(f_extra)
    if seqs:
        # <----------- MAX ------------>
        # [1][2][3][4][5][6][7][8][9][X] # MAX_PROT
        #                               xxxxxxxxxxxx    # untranslated_overhang
        # .....................[E][P][A]............
        # .....................[E][P][A].............
        # .....................[E][P][A]..............
        #                               xxxxxxxxxxxxxx  # max_untranslated_overhang

        with open(f_seq_nuc, 'w', encoding='utf-8') as f, open(
            f_seq_aa, 'w', encoding='utf-8'
        ) as f_aa:
            if not MAX:
                MAX = 12
            MAX = max(
                12,
                MAX,
                (stats['MAX_LEN'] if stats else max(map(len, seqs)))
                - max_untranslated_overhang,
            )   # ensure at least 4 aa columns
            MAX_PROT = MAX // 3
            # print(MAX, MAX_PROT)
            f.write(
                f'#{" "*(MAX-10)}[E][P][A]{" "*max_untranslated_overhang}\n'
            )
            f_aa.write(f'#{" "*(MAX_PROT-4)}EPA\n')
            for seq in seqs:
                f.write(f'{seq: >{MAX+max_untranslated_overhang}}\n')
                f_aa.write(
                    f'{seq2aa(seq[:3]).lower()+seq2aa(seq[3:-max_untranslated_overhang]): >{MAX_PROT}}\n'
                )   # 14 nuc after A-site


def export_all(results_all, outdir=None):
    """Exports files for all results using :func:`export_data`"""
    for filename, (seqs, stats) in results_all.items():
        print(filename)
        print(simple_graph(stats['lengths'], start=0, wrap=81))
        export_data(filename, seqs=seqs, stats=stats, outdir=outdir, MAX=None)


@wraps(parse_filter_fastq)   # FIXME add parameters of export_data
def parse(
    filename,
    *,
    seqs=None,
    stats=None,
    outdir=None,
    MAX=None,
    untranslated_overhang=12,
    **kwargs,
):
    """Wrapper to combine parse_filter_fastq and export_data"""
    seqs, stats, _ = parse_filter_fastq(filename, **kwargs)
    export_data(
        filename,
        seqs=seqs,
        stats=stats,
        outdir=outdir,
        MAX=MAX,
        untranslated_overhang=untranslated_overhang,
    )


def format_sequences(
    filename, codons=False, aa=False, repeat_header=False, out=None, limit=None
):
    """
    Formats a nucleotide inverse toeprint file in a custom human-readable format

    Parameters
    ----------

    filename : str or list
        Name of the nucleotide inverse toeprint file to use as input.
        If a list of filenames or a directory is passed, apply the function to all the nucleotide inverse toeprint files.
    codons : bool
        If True, splits the coding sequence into codons
    aa: bool
        If True, interleave the codons below each nucleotide sequence
        Adds the length of the peptide after the amino-acids.
    repeat_header : bool or int:
        If given an integer, repeats the header every <repeat_sequence> reads.
    out : str or Path
        If defined, write the output to this file. If None, write to stdout.
    limit: int or None
        Limit the number of reads to process

    Examples
    --------
    Display the first 5 inverse toeprints of the "nnn15_noa1.nuc.itp.txt" file:

     >>> format_sequences('nnn15_noa1.nuc.itp.txt', limit=5)
     #                    [E][P][A]
                          ATGGGACGC cccgcagtatct
           ATGAGTTACAAAGGCAACTCGGAA caggtagcatatc
                          ATGGAAGAG gcccatgccattcc
                             ATGAAT cgaaacatgttt
     ATGACTATGTTTCTTGGACACACATAAGGG aactagttaggg

    Display the first 10 inverse toeprints of the "nnn15_noa1.nuc.itp.txt" file,
    group the coding sequence by codons and display the translation.
    Repeat the header every 5 reads.:

    .. code-block:: python

         >>> format_sequences('nnn15_noa1.nuc.itp.txt', limit=5,
         ...                  codons=True, aa=True,
         ...                  repeat_header=5)
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

         #                           [E] [P] [A]
                                     ATG CTA TAA taggtcaagcacca
                                      M   L   *  3
                             ATG ACC AAT CCG TAG gactaacgccacat
                              M   T   N   P   *  5
                     ATG TAG CCG GGC AAG GAG ATC cgcacctcgcgc
                      M   *   P   G   K   E   I  7
                                         ATG TAA ctatacgacgtcg
                                          M   *  2
                                         ATG TAA acacgccttgtcgt
                                          M   *  2

    Export the output to a file

    >>> format_sequences('nnn15_noa1.nuc.itp.txt', out)
    """

    SUFFIX = f'.nuc.{ITP_FILE_SUFFIX}.txt'

    if out is not None:
        out = Path(out)

    # if out is a list of files or a directory, apply on all the files
    if (is_list := isinstance(filename, list)) or (
        dir_ := Path(filename)
    ).is_dir():
        if not is_list:
            filename = sorted(list(dir_.glob(f'*{SUFFIX}')))
        if out is not None:
            if (out := Path(out)).exists():
                if not out.is_dir():
                    raise ValueError(
                        'If several files or a directory are passed as input, "out" should be a directory.'
                    )
            else:
                out.mkdir(parents=True)
        for f in filename:
            if out is None:
                # if no output dir, display the input file name before each output
                print(f'{f}')
            format_sequences(
                f,
                codons=codons,
                aa=aa,
                repeat_header=repeat_header,
                limit=limit,
                out=out,
            )
        return   # delegate processing to a new function call, stop here

    if isinstance(out, Path) and out.is_dir():
        # If out is a directory, create a filename based on the original name
        out /= Path(filename).name.removesuffix(SUFFIX) + f'.formatted{SUFFIX}'

    with open(filename) as f, nullcontext(sys.stdout) if out is None else open(
        out, 'w'
    ) as f_out:
        first = f.readline()
        leading = len(first.rstrip())
        f.seek(0)

        for i, s in enumerate(f):
            s1 = s[:leading]
            pep_len = len(s1.strip()) // 3
            if codons:
                s1 = ' '.join([s[i : i + 3] for i in range(0, leading, 3)])
            s2 = f'{s1} {s[leading:].lower()}'
            f_out.write(s2)

            if aa and not s1.startswith('#'):
                sa = (' ' if codons else '').join(
                    [
                        f' {codon_table.get(s[i:i+3], " ")} '
                        for i in range(0, leading, 3)
                    ]
                )
                f_out.write(f'{sa} {pep_len}\n')

            if i == 0:
                header = '\n' + s2
            elif repeat_header and i < limit and not i % repeat_header:
                f_out.write(header)

            if limit and i >= limit:
                break


def main(out=False):
    """Wrapper for the command line interface"""
    import argparse
    import datetime

    def ArgRange(string):
        import re
        from argparse import ArgumentTypeError

        if string.isdigit():
            return (int(string),) * 2

        m = re.match(r'^(?:(\d+)?-)?(\d+)?$', string)
        if not m:
            raise ArgumentTypeError(
                f'invalid format "{string}", use 1-10 or 1- or -10 or 10'
            )
        return tuple(int(x) if x else 0 for x in m.groups())

    parser = argparse.ArgumentParser(
        description='parse fastq files and extract peptides sequence for inverse toe-printing',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('files', help='input files', nargs='+')

    parser.add_argument(
        '-o',
        '--outdir',
        help='directory for the output files',
        dest='outdir',
        default='.',
    )

    parser.add_argument(
        '-a1',
        '--left-adaptor',
        help='Sequence (5ʹ→3ʹ) of the left adaptor',
        dest='a1',
        default=DEFAULTS['a1'],
    )
    parser.add_argument(
        '-a2',
        '--right-adaptor',
        help='Sequence (5ʹ→3ʹ) of the right adaptor',
        dest='a2',
        default=DEFAULTS['a2'],
    )

    parser.add_argument(
        '-s',
        '--peptide-size',
        help='Allowed length range of peptide considered (min-max)',
        dest='range',
        default='1-10',
        type=ArgRange,
    )
    parser.add_argument(
        '-q',
        '--quality',
        help='Threshold for the PHRED quality cutoff',
        dest='quality',
        default=30,
    )

    parser.add_argument(
        '--limit',
        help='Max number of sequences to process (useful for quick tests)',
        dest='limit',
        type=int,
        default=None,
    )

    args = parser.parse_args()

    args_dic = vars(args)
    from pprint import pprint

    pprint(args_dic)
    kwargs = args_dic.copy()

    kwargs['min_seq_len'] = args_dic['range'][0] * 3
    max_peptide = args_dic['range'][1]
    kwargs['max_seq_len'] = max_peptide * 3 if max_peptide else None

    del kwargs['files']
    del kwargs['outdir']
    print(datetime.datetime.now().strftime('Start time: %Y-%m-%d %H:%M'))
    results_all = parse_all(
        args_dic['files'], outdir=args_dic['outdir'], save=True, **kwargs
    )
    print(datetime.datetime.now().strftime('End time: %Y-%m-%d %H:%M'))
    if out:
        return results_all


if __name__ == '__main__':
    main()
