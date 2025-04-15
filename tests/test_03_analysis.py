import os
from pathlib import Path

import pytest
from itpseq import DataSet, Replicate, Sample, processing


class TestFunctions:
    def test_read_itp_file_as_series_nuc(self, data_dir):
        ser = processing.read_itp_file_as_series(
            data_dir / 'tcx_small_test' / 'nnn15_noa1.nuc.itp.txt',
            how='aax',
        )
        assert ser.head(10).tolist() == [
            '      ATGAGTTACAAAGGCAACTCGGAACAGGTAGCATATC ',
            '                     ATGGAAGAGGCCCATGCCATTCC',
            '                     ATGCTATAATAGGTCAAGCACCA',
            '               ATGACCAATCCGTAGGACTAACGCCACAT',
            '                        ATGTAACTATACGACGTCG ',
            '                        ATGTAAACACGCCTTGTCGT',
            '                        ATGACACCTACAGGGCCCAT',
            '                  ATGACAAAATAGCAGCGCAACGGGA ',
            '                     ATGAACCGGTAACAAGTCCACCT',
            '                        ATGTCGGGATCGATGCAGA ',
        ]

    def test_read_itp_file_as_series_aax(self, data_dir):
        ser = processing.read_itp_file_as_series(
            data_dir / 'tcx_small_test' / 'nnn15_noa1.aa.itp.txt',
            how='aax',
        )
        assert ser.head(10).tolist() == [
            '       mGR',
            '  mSYKGNSE',
            '       mEE',
            '        mN',
            '       mL*',
            '     mTNP*',
            '        m*',
            '        m*',
            '        mT',
            '     mIRD*',
        ]

    def test_read_itp_file_as_series_aa(self, data_dir):
        ser = processing.read_itp_file_as_series(
            data_dir / 'tcx_small_test' / 'nnn15_noa1.aa.itp.txt',
            how='aa',
        )
        assert ser.head(10).tolist() == [
            '       mGR',
            '  mSYKGNSE',
            '       mEE',
            '        mN',
            'mTMFLGHT*G',
            '       mL*',
            '     mTNP*',
            '   m*PGKEI',
            '        m*',
            '        m*',
        ]

    def test_read_itp_file_as_series_limit(self, data_dir):
        ser = processing.read_itp_file_as_series(
            data_dir / 'tcx_small_test' / 'nnn15_noa1.aa.itp.txt',
            limit=3,
        )
        assert ser.head(10).tolist() == [
            '       mGR',
            '  mSYKGNSE',
            '       mEE',
        ]

    def test_read_itp_file_as_series_sample(self, data_dir):
        ser = processing.read_itp_file_as_series(
            data_dir / 'tcx_small_test' / 'nnn15_noa1.aa.itp.txt',
            sample=3,
        )
        assert len(ser) == 3

    @pytest.mark.parametrize(
        'inpt, params, expected',
        [
            ('-2', {}, -4),
            ('E', {}, -3),
            ('P', {}, -2),
            (0, {}, -2),
            ('A', {}, -1),
            (0, {'how': 'nuc'}, -20),
        ],
    )
    def test_code2pos(self, inpt, params, expected):
        assert processing.code2pos(inpt, **params) == expected

    @pytest.mark.parametrize(
        'inpt, params',
        [
            ('X', {}),
            ('X', {'how': 'nuc'}),
        ],
    )
    def test_code2pos_err(self, inpt, params):
        with pytest.raises(ValueError):
            processing.code2pos(inpt, **params)

    @pytest.mark.parametrize(
        'inpt, params, expected',
        [
            ('P', {}, [-2]),
            ('E,A', {}, [-3, -1]),
            ('-3,E:A', {}, [-5, slice(-3, None, None)]),
            ('E:P', {}, [slice(-3, -1, None)]),
            ('E:A', {}, [slice(-3, None, None)]),
            ('-1:1', {}, [slice(-3, None, None)]),
            ('-2:2', {'how': 'nuc'}, [slice(-22, -17, None)]),
            ('-2,0', {'how': 'nuc'}, [-22, -20]),
        ],
    )
    def test_ranges(self, inpt, params, expected):
        assert processing.ranges(inpt, **params) == expected

    @pytest.mark.parametrize(
        'inpt, expected',
        [
            (-2, '-6:-4'),
            ('-2', '-6:-4'),
            ('E', '-3:-1'),
            ('P', '0:2'),
            ('A', '3:5'),
            ('-2:P', '-6:2'),
            ('-2:E,A', '-6:-1,3:5'),
        ],
    )
    def test_aa2nuc_pos(self, inpt, expected):
        assert processing.aa2nuc_pos(inpt) == expected


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


@pytest.mark.xdist_group(name='sample_analysis')
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

    def test_get_counts_ratio_EP(self):
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

    def test_get_counts_ratio_nuc(self):
        self.tcx_data._clear_cache(force=True)
        d = (
            self.tcx_data['nnn15.tcx']
            .get_counts_ratio(pos='-3:2', how='nuc')
            .sum()
            .astype(int)
            .to_dict()
        )
        assert d == {
            'nnn15.noa.1': 40805,
            'nnn15.noa.2': 34527,
            'nnn15.noa.3': 38217,
            'nnn15.tcx.1': 37999,
            'nnn15.tcx.2': 27134,
            'nnn15.tcx.3': 33692,
            'nnn15.noa': 1009599,
            'nnn15.tcx': 1013385,
            'ratio': 4876,
        }

    def test_get_counts_ratio_pos(self):
        self.tcx_data._clear_cache(force=True)
        df = self.tcx_data['nnn15.tcx'].get_counts_ratio_pos()
        assert set(df.index) == {'-2', 'E', 'P', 'A'}
        assert set(df.columns) == set('*ACDEFGHIKLMNPQRSTVWYm')
        assert df.notna().values.sum() == 85

    def test_get_counts_ratio_pos_nuc(self):
        self.tcx_data._clear_cache(force=True)
        df = self.tcx_data['nnn15.tcx'].get_counts_ratio_pos(how='nuc')
        assert set(df.index) == set(range(-6, 6))
        assert set(df.columns) == set('ACGT')
        assert df.notna().values.sum() == 48

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

        df = self.tcx_data['nnn15.tcx'].DE('E:P', join=True)
        assert list(df.columns) == [
            'nnn15.noa.1',
            'nnn15.noa.2',
            'nnn15.noa.3',
            'nnn15.tcx.1',
            'nnn15.tcx.2',
            'nnn15.tcx.3',
            'nnn15.noa',
            'nnn15.tcx',
            'ratio',
            'baseMean',
            'log2FoldChange',
            'lfcSE',
            'stat',
            'pvalue',
            'padj',
            'log10pvalue',
            'log10padj',
        ]
