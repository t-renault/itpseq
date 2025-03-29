import os
from pathlib import Path

import pytest
from itpseq import DataSet, Replicate, Sample


class TestDataSet:
    @pytest.fixture(scope='class', autouse=True)
    def setup_class(self, request, data_dir, tmp_outdir):
        request.cls.tcx_data = DataSet(
            data_dir / 'tcx_small_test', result_path=tmp_outdir
        )

    def test_infos(self):
        self.tcx_data._clear_cache(force=True)
        df = self.tcx_data.infos()
        assert set(df.index) == {
            'nnn15.noa.1',
            'nnn15.noa.2',
            'nnn15.noa.3',
            'nnn15.tcx.1',
            'nnn15.tcx.2',
            'nnn15.tcx.3',
        }
        assert list(df.columns) == [
            'total_sequences',
            'noadaptor',
            'contaminant',
            'lowqual',
            'tooshort',
            'toolong',
            'extra0',
            'extra1',
            'extra2',
            'MAX_LEN',
        ]

    def test_toeprint_df(self):
        self.tcx_data._clear_cache(force=True)
        df = self.tcx_data.toeprint_df
        assert set(df.columns) == {
            'nnn15.noa.1',
            'nnn15.noa.2',
            'nnn15.noa.3',
            'nnn15.tcx.2',
            'nnn15.tcx.3',
            'nnn15.tcx.1',
        }
        # needs to be updated if the dataset changes
        assert df.shape == (98, 6)
        assert df.sum().sum() == 551145


# noinspection PyComparisonWithNone
class TestSample:
    @pytest.fixture(scope='class', autouse=True)
    def setup_class(self, request, data_dir):
        request.cls.tcx_data = DataSet(data_dir / 'tcx_small_test')

    def test_copy(self):
        self.tcx_data._clear_cache(force=True)
        new = self.tcx_data['nnn15.tcx'].copy(name='new', reference=None)
        assert new.name == 'new'
        assert new.reference is None

    def test_infos(self):
        self.tcx_data._clear_cache(force=True)
        df = self.tcx_data['nnn15.tcx'].infos()
        assert set(df.index) == {
            'nnn15.tcx.1',
            'nnn15.tcx.2',
            'nnn15.tcx.3',
        }
        assert list(df.columns) == [
            'total_sequences',
            'noadaptor',
            'contaminant',
            'lowqual',
            'tooshort',
            'toolong',
            'extra0',
            'extra1',
            'extra2',
            'MAX_LEN',
        ]

    def test_itp_len(self):
        self.tcx_data._clear_cache(force=True)
        df = self.tcx_data['nnn15.tcx'].itp_len
        assert list(df.columns) == ['length', 'replicate', 'count', 'sample']
        assert len(df) > 200

    def test_get_counts(self):
        self.tcx_data._clear_cache(force=True)
        df = self.tcx_data['nnn15.tcx'].get_counts()
        assert df.shape == (23, 30)
        assert tuple(map(set, df.columns.levels)) == (
            {'nnn15.tcx.1', 'nnn15.tcx.2', 'nnn15.tcx.3'},
            {-8, -7, -6, -5, -4, -3, -2, -1, 0, 1},
        )
        # fmt: off
        assert set(df.index) == {
            ' ', '*', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'm',
        }
        # fmt: on
        # needs to be updated if the dataset changes
        assert df.sum().sum() == 849170

    def test_get_counts_E(self):
        self.tcx_data._clear_cache(force=True)
        df = self.tcx_data['nnn15.tcx'].get_counts('E')
        assert df.shape == (22, 3)
        assert set(df.columns) == {'nnn15.tcx.1', 'nnn15.tcx.2', 'nnn15.tcx.3'}
        # fmt: off
        assert set(df.index) == {
            ' ', 'm', 'T', 'P', 'S', 'R', 'L', 'A', 'Q', 'K', 'H',
            'V', 'I', 'N', 'G', 'E', 'Y', 'D', 'F', 'C', 'W', 'M',
        }
        # fmt: on
        # needs to be updated if the dataset changes
        assert df.sum().sum() == 84917

    def test_get_counts_EA(self):
        self.tcx_data._clear_cache(force=True)
        df = self.tcx_data['nnn15.tcx'].get_counts('E:A')
        assert df.shape == (7640, 3)
        assert set(df.columns) == {'nnn15.tcx.1', 'nnn15.tcx.2', 'nnn15.tcx.3'}
        # needs to be updated if the dataset changes
        assert df.sum().sum() == 84917

    def test_get_counts_EP(self):
        self.tcx_data._clear_cache(force=True)
        df = self.tcx_data['nnn15.tcx'].get_counts_ratio('E:P')
        assert df.shape == (421, 9)
        assert set(df.columns) == {
            'nnn15.noa.1',
            'nnn15.noa.2',
            'nnn15.noa.3',
            'nnn15.tcx.1',
            'nnn15.tcx.2',
            'nnn15.tcx.3',
            'nnn15.noa',
            'nnn15.tcx',
            'ratio',
        }

    def test_get_DE_EP(self):
        self.tcx_data._clear_cache(force=True)
        df = self.tcx_data['nnn15.tcx'].DE('E:P')
        assert df.shape == (420, 8)
        assert list(df.columns) == [
            'baseMean',
            'log2FoldChange',
            'lfcSE',
            'stat',
            'pvalue',
            'padj',
            'log10pvalue',
            'log10padj',
        ]
