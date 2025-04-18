import filecmp
import sys
from collections import Counter
from itertools import product
from unittest.mock import patch

import itpseq.parsing
import pytest


class TestParse:
    def test_parse_all(self, data_dir, tmp_outdir):
        path = data_dir / 'tcx_small_test'
        d = itpseq.parsing.parse_all(
            pattern=str(path / '*.fastq'), outdir=tmp_outdir, save=True
        )
        # check all files were processed
        assert set(d) == {
            'nnn15_noa1.assembled.fastq',
            'nnn15_noa2.assembled.fastq',
            'nnn15_noa3.assembled.fastq',
            'nnn15_tcx1.assembled.fastq',
            'nnn15_tcx2.assembled.fastq',
            'nnn15_tcx3.assembled.fastq',
        }
        files = [
            ''.join(p)
            for p in product(
                (k.removesuffix('.assembled.fastq') for k in d),
                ('.aa.itp.txt', '.nuc.itp.txt', '.itp.log', '.itp.json'),
            )
        ]
        # check the content of the output files is the same
        ok, nok, err = filecmp.cmpfiles(path, tmp_outdir, files, shallow=False)
        assert nok == [] and err == []

    @pytest.mark.parametrize(
        'ext', ['', '.gz', '.bz2', '.xz', '.lzma', '.zip']
    )
    def test_cli_parse(self, ext, data_dir, tmp_outdir):
        """Extra test due to an issue to correctly report test coverage with multiprocessing"""
        path = data_dir / 'tcx_small_test'
        itpseq.parsing.parse(
            path / f'nnn15_noa1.assembled.fastq{ext}', outdir=tmp_outdir
        )
        assert filecmp.cmp(
            path / 'nnn15_noa1.nuc.itp.txt',
            tmp_outdir / 'nnn15_noa1.nuc.itp.txt',
            shallow=True,
        )

    @pytest.mark.parametrize(
        'file_prefix, params',
        [
            ('nnn15_noa1', ()),
            ('nnn15_noa2', ('-s', '1-10')),
        ],
    )
    def test_legacy_parse_cli(self, data_dir, tmp_outdir, file_prefix, params):
        path = data_dir / 'tcx_small_test'
        input_file = str(path / f'{file_prefix}.assembled.fastq')
        expected_file_patterns = [
            f'{file_prefix}{suffix}'
            for suffix in (
                '.aa.itp.txt',
                '.nuc.itp.txt',
                '.itp.log',
                '.itp.json',
            )
        ]
        test_args = ['itpseq', input_file, '-o', str(tmp_outdir), *params]
        with patch.object(sys, 'argv', test_args):
            result = itpseq.parsing.main(out=True)

        ok, nok, err = filecmp.cmpfiles(
            path, tmp_outdir, expected_file_patterns, shallow=False
        )

        assert (
            nok == [] and err == []
        ), f"Files don't match: nok={nok}, err={err}"
        assert result is not None

    def test_parse_limit_max_seq_len(self, data_dir, tmp_outdir):
        path = data_dir / 'tcx_small_test'
        itpseq.parsing.parse(
            path / f'nnn15_noa1.assembled.fastq',
            outdir=tmp_outdir,
            limit=20,
            max_seq_len=40,
        )
        with open(tmp_outdir / 'nnn15_noa1.nuc.itp.txt') as f:
            assert f.readlines() == [
                '#                             [E][P][A]              \n',
                '                              ATGGGACGCCCCGCAGTATCT  \n',
                '               ATGAGTTACAAAGGCAACTCGGAACAGGTAGCATATC \n',
                '                              ATGGAAGAGGCCCATGCCATTCC\n',
                'ATGTCTGCGAGATGATCACGGATGATCGTAGACGAACTACACCCGACTGCG  \n',
                '                                 ATGAATCGAAACATGTTT  \n',
                '         ATGACTATGTTTCTTGGACACACATAAGGGAACTAGTTAGGG  \n',
            ]

    def test_format_sequences(self, data_dir, tmp_outdir):
        path = data_dir / 'tcx_small_test'
        out = tmp_outdir / 'out.txt'
        itpseq.parsing.format_sequences(
            path / 'nnn15_noa1.nuc.itp.txt',
            codons=True,
            aa=True,
            repeat_header=2,
            out=out,
            limit=3,
        )
        with open(out) as f:
            assert f.readlines() == [
                '#                           [E] [P] [A]               \n',
                '                            ATG GGA CGC cccgcagtatct  \n',
                '                             M   G   R  3\n',
                '        ATG AGT TAC AAA GGC AAC TCG GAA caggtagcatatc \n',
                '         M   S   Y   K   G   N   S   E  8\n',
                '\n',
                '#                           [E] [P] [A]               \n',
                '                            ATG GAA GAG gcccatgccattcc\n',
                '                             M   E   E  3\n',
            ]

    def test_format_sequences_dir(self, capsys, monkeypatch, data_dir):
        # temporarily change directory for this test
        monkeypatch.chdir(data_dir / 'tcx_small_test')
        itpseq.parsing.format_sequences('.', limit=1)
        # capture stdout
        captured = capsys.readouterr()
        assert captured.out.splitlines() == [
            'nnn15_noa1.nuc.itp.txt',
            '#                    [E][P][A]               ',
            '                     ATGGGACGC cccgcagtatct  ',
            'nnn15_noa2.nuc.itp.txt',
            '#                    [E][P][A]               ',
            '                        ATGTAC tcacatccccgatc',
            'nnn15_noa3.nuc.itp.txt',
            '#                    [E][P][A]               ',
            '         ATGGGCCTATACTGTTATGAA ctacccatagcagt',
            'nnn15_tcx1.nuc.itp.txt',
            '#                    [E][P][A]               ',
            '                        ATGTCA cgtcaaacccaagt',
            'nnn15_tcx2.nuc.itp.txt',
            '#                    [E][P][A]               ',
            '               ATGCACATAACCTAA tgaatcctaactc ',
            'nnn15_tcx3.nuc.itp.txt',
            '#                    [E][P][A]               ',
            '                     ATGAATACG tattgtcaaccct ',
        ]

    @pytest.mark.parametrize(
        'inpt, params, expected',
        [
            (Counter(), {}, ''),
            (
                Counter({0: 100, 10: 50, 20: 200}),
                {},
                '▄_________▂_________█\n╹         ╹\n0        10',
            ),
            (
                Counter({0: 100, 10: 50, 20: 200}),
                {'ticks': 5},
                '▄_________▂_________█\n╹    ╹    ╹    ╹\n0    5   10   15',
            ),
            (
                Counter({0: 100, 10: 50, 20: 200}),
                {'wrap': 10, 'ticks': 5},
                '▄_________\n╹    ╹\n0    5\n▂_________\n╹    ╹\n10   15',
            ),
            (
                Counter({0: 100, 10: 50, 20: 200}),
                {'vmax': 100},
                '█_________▄_________▒\n╹         ╹\n0        10',
            ),
        ],
    )
    def test_simple_graph(self, inpt, params, expected):
        assert itpseq.parsing.simple_graph(inpt, **params) == expected
