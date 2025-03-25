import pytest
import os
from pathlib import Path

from itpseq import DataSet, Sample, Replicate

data_dir = Path(__file__).parent / 'test_data/'


class TestDataSet:
    def test_DataSet_from_path(self, tmp_outdir):
        data = DataSet(data_dir / 'tcx_small_test/', result_path=tmp_outdir)
        assert set(data.samples) == {'nnn15.noa', 'nnn15.tcx'}
        assert set(data.samples_with_ref) == {'nnn15.tcx'}
        assert set(data.replicates) == {
            'nnn15.noa.1',
            'nnn15.noa.2',
            'nnn15.noa.3',
            'nnn15.tcx.1',
            'nnn15.tcx.2',
            'nnn15.tcx.3',
        }
        assert data.keys == ['lib_type', 'sample']
        assert data.ref_labels == (('sample', 'noa'),)
        assert data['nnn15.tcx'].labels == {'lib_type': 'nnn15', 'sample': 'tcx'}
        assert data['nnn15.noa'].labels == {'lib_type': 'nnn15', 'sample': 'noa'}
        assert repr(data) == """DataSet(data_path='"""+str(data_dir)+"""/tcx_small_test',
        file_pattern='(?P<lib_type>[^_]+)_(?P<sample>[^_\\\\d]+)(?P<replicate>\\\\d+)',
        samples=[Sample(nnn15.noa:[1, 2, 3]),
                 Sample(nnn15.tcx:[1, 2, 3], ref: nnn15.noa)],
        )"""

    def test_DataSet_from_path_keys(self, tmp_outdir):
        data = DataSet(
            data_dir / 'tcx_small_test/',
            keys=['sample'],
            result_path=tmp_outdir,
        )
        assert set(data.samples) == {'noa', 'tcx'}
        assert set(data.samples_with_ref) == {'tcx'}
        assert set(data.replicates) == {
            'noa.1',
            'noa.2',
            'noa.3',
            'tcx.1',
            'tcx.2',
            'tcx.3',
        }
        assert data.keys == ['sample']
        assert data.ref_labels == (('sample', 'noa'),)

    def test_DataSet_from_path_pattern(self, tmp_outdir):
        data = DataSet(
            data_dir / 'loading_concentrations/',
            file_pattern=r'(?P<sample>[^_]+)(?P<replicate>\d+)(_(?P<concentration>\d+µM))?',
            result_path=tmp_outdir,
        )
        assert set(data.samples) == {
            'drugA.10µM',
            'drugA.20µM',
            'drugA.30µM',
            'drugB.10µM',
            'drugB.20µM',
            'drugB.30µM',
            'noa',
        }
        assert set(data.samples_with_ref) == {
            'drugA.10µM',
            'drugA.20µM',
            'drugA.30µM',
            'drugB.10µM',
            'drugB.20µM',
            'drugB.30µM',
        }

    def test_DataSet_from_dictionary(self, tmp_outdir):
        data = DataSet(
            {
                'tcx': [
                    {'file_prefix': 'tcx_small_test/nnn15_tcx1'},
                    {'file_prefix': 'tcx_small_test/nnn15_tcx2'},
                    {'file_prefix': 'tcx_small_test/nnn15_tcx3'},
                ],
                'noa': [
                    {'file_prefix': 'tcx_small_test/nnn15_noa1'},
                    {
                        'file_prefix': 'tcx_small_test/nnn15_noa2',
                        'replicate': 'custom_name',
                    },
                    {'file_prefix': 'tcx_small_test/nnn15_noa3'},
                ],
            },
            ref_mapping={'tcx': 'noa'},
            result_path=tmp_outdir,
        )
        assert set(data.samples) == {'noa', 'tcx'}
        assert set(data.replicates) == {
            'tcx.rep1',
            'tcx.rep2',
            'tcx.rep3',
            'noa.rep1',
            'noa.custom_name',
            'noa.rep3',
        }
