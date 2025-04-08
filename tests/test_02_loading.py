import os
from pathlib import Path

import pytest
from itpseq import DataSet, Replicate, Sample

data_dir = Path(__file__).parent / 'test_data'


class TestDataSet:
    def test_dataset_from_path(self, tmp_outdir):
        data = DataSet(data_dir / 'tcx_small_test', result_path=tmp_outdir)
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
        # fmt: off
        assert data['nnn15.tcx'].labels == {'lib_type': 'nnn15', 'sample': 'tcx'}
        assert data['nnn15.noa'].labels == {'lib_type': 'nnn15', 'sample': 'noa'}

        # testing __repr__ and __gt__
        assert repr(data) == """DataSet(data_path='"""+str(data_dir/'tcx_small_test')+"""',
        file_pattern='(?P<lib_type>[^_]+)_(?P<sample>[^_\\\\d]+)(?P<replicate>\\\\d+)',
        samples=[Sample(nnn15.noa:[1, 2, 3]),
                 Sample(nnn15.tcx:[1, 2, 3], ref: nnn15.noa)],
        )"""
        # fmt: on
        assert repr(data['nnn15.noa']) == 'Sample(nnn15.noa:[1, 2, 3])'
        assert (
            repr(data['nnn15.tcx'])
            == 'Sample(nnn15.tcx:[1, 2, 3], ref: nnn15.noa)'
        )
        assert repr(data.replicates['nnn15.noa.1']) == 'Replicate(nnn15.noa.1)'
        assert data['nnn15.noa'] < data['nnn15.tcx']

        assert data.replicates['nnn15.noa.1'].sample_name == 'nnn15.noa'
        assert data.replicates['nnn15.noa.1'].filename == str(
            data_dir / 'tcx_small_test' / 'nnn15_noa1.nuc.itp.txt'
        )

    def test_dataset_from_path_keys(self, tmp_outdir):
        data = DataSet(
            data_dir / 'tcx_small_test',
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

    def test_dataset_from_path_pattern(self, tmp_outdir):
        data = DataSet(
            data_dir / 'loading_concentrations',
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

    def test_dataset_from_dictionary(self, tmp_outdir):
        data = DataSet(
            {
                'tcx': [
                    {'file_prefix': Path('tcx_small_test') / 'nnn15_tcx1'},
                    {'file_prefix': Path('tcx_small_test') / 'nnn15_tcx2'},
                    {'file_prefix': Path('tcx_small_test') / 'nnn15_tcx3'},
                ],
                'noa': [
                    {'file_prefix': Path('tcx_small_test') / 'nnn15_noa1'},
                    {
                        'file_prefix': Path('tcx_small_test') / 'nnn15_noa2',
                        'replicate': 'custom_name',
                    },
                    {'file_prefix': Path('tcx_small_test') / 'nnn15_noa3'},
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

    def test_dataset_from_json(self, tmp_outdir):
        import json

        with open(data_dir / 'tcx_small_test' / 'samples.json') as f:
            data = DataSet(json.load(f))
        assert set(data.samples) == {'noa', 'tcx'}
        assert set(data.replicates) == {
            'tcx.rep1',
            'tcx.rep2',
            'tcx.rep3',
            'noa.rep1',
            'noa.custom_name',
            'noa.rep3',
        }


class TestSample:
    def test_sample_empty(self):
        s = Sample()
        assert len(s.name) == 36       # check UUID
        assert s.name.count('-') == 4  #
        assert s.replicates == []

    def test_sample_from_list_of_dicts(self):
        noa = Sample(
            [
                {'file_prefix': 'nnn15_noa1'},
                {'file_prefix': 'nnn15_noa2'},
                {'file_prefix': 'nnn15_noa3'},
            ],
            name='noa',
        )
        tcx = Sample(
            [
                {'file_prefix': 'nnn15_tcx1'},
                {'file_prefix': 'nnn15_tcx2'},
                {'file_prefix': 'nnn15_tcx3'},
            ],
            name='tcx',
            reference=noa,
        )
        assert repr(tcx) == 'Sample(tcx:[rep1, rep2, rep3], ref: noa)'
        assert (
            repr(noa.replicates)
            == '[Replicate(noa.rep1), Replicate(noa.rep2), Replicate(noa.rep3)]'
        )
        assert tcx > noa
        assert tcx.copy(name='abc', reference=noa) > noa
        assert tcx.copy(name='abc') < noa
        assert tcx.copy(name='xyz') > noa


class TestReplicate:
    def test_replicate_gt(self):
        assert Replicate(replicate=1) < Replicate(replicate=2)
        assert Replicate(replicate=2) > Replicate(replicate=1)
        assert Replicate(replicate='2') < Replicate(replicate='10')
        assert Replicate(replicate=2) < Replicate(replicate='10')
        assert Replicate(replicate='A') < Replicate(replicate='Z')
