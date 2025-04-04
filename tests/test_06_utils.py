import pytest
from itpseq import utils


class TestUtils:
    def test_itp_schematic(self):
        params = [
            {},
            {'codons': 2},
            {'coding_sequence': 'ATTACTCCGAGCGAGCAA'},
        ]
        checks = [
            'NNNNNNNNNNNNNNNNNNNNNNNNNNNnnnnnnnnnnnn--mXXXXXXXX',
            'NNNNNNnnnnnnnnnnnn--mX',
            'ATTACTCCGAGCGAGCAAnnnnnnnnnnnn--ITPSEQ',
        ]
        for p, c in zip(params, checks):
            ax = utils.itp_schematic(**p)
            # check that E P A overhang patches are drawn
            assert len(ax.patches) == 4
            # checks the texts (patches, nucleotides, amino-acids)
            texts = [
                c.get_text()
                for c in ax.get_children()
                if hasattr(c, 'get_text')
            ]
            assert texts[:4] == ['protected overhang', 'E', 'P', 'A']
            assert ''.join(texts[4:]) == c

    def test_dict_to_tuple(self):
        tests = [
            # inpt, params, out
            ({}, {}, ()),
            ({'A': 'a'}, {}, (('A', 'a'),)),
            (
                {'A': 'a', 'B': 'b', 'C': 'c'},
                {'ignore': 'A'},
                (('B', 'b'), ('C', 'c')),
            ),
            ({'A': 'a', 'B': 'b', 'C': 'c'}, {'keep': 'B'}, (('B', 'b'),)),
            (
                {'A': 'a', 'B': 'b', 'C': 'c'},
                {'keep': ['B', 'C']},
                (('B', 'b'), ('C', 'c')),
            ),
        ]
        for inpt, params, out in tests:
            assert utils.dict_to_tuple(inpt, **params) == out

    def test_dedup_names(self):
        tests = [
            # inpt, params, out
            ([], {}, []),
            (['x', 'y', 'x', 'x'], {}, ['x', 'y', 'x.1', 'x.2']),
            (['x', 'y', 'x', 'x'], {'sep': '-'}, ['x', 'y', 'x-1', 'x-2']),
            (['x', 'y', 'x.1', 'x'], {}, ['x', 'y', 'x.1', 'x.1.1']),
        ]
        for inpt, params, out in tests:
            assert utils.dedup_names(inpt, **params) == out
