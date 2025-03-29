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
