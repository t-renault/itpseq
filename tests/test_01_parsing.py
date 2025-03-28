import filecmp
from itertools import product
from pathlib import Path

import pytest

import itpseq.parsing


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
        ok, nok, err = filecmp.cmpfiles(path, tmp_outdir, files)
        assert nok == [] and err == []
