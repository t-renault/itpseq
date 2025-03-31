import filecmp
from itertools import product
from pathlib import Path

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

    def test_cli_parse(self, data_dir, tmp_outdir):
        """Extra test due to an issue to correctly report test coverage with multiprocessing"""
        path = data_dir / 'tcx_small_test'
        itpseq.parsing.parse(path / 'nnn15_noa1.assembled.fastq')

    def test_format_sequences(self, data_dir, tmp_outdir):
        path = data_dir / 'tcx_small_test'
        out = Path(tmp_outdir) / 'out.txt'
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
