import pytest
from itpseq import utils


class TestUtils:
    @pytest.mark.parametrize(
        'params, expected',
        [
            ({}, 'NNNNNNNNNNNNNNNNNNNNNNNNNNNnnnnnnnnnnnn--mXXXXXXXX'),
            ({'codons': 2}, 'NNNNNNnnnnnnnnnnnn--mX'),
            (
                {'coding_sequence': 'ATTACTCCGAGCGAGCAA'},
                'ATTACTCCGAGCGAGCAAnnnnnnnnnnnn--ITPSEQ',
            ),
            (
                {'coding_sequence': 'ATG', 'overhang_sequence': 'test-'},
                'ATGtest-M',
            ),
            (
                {'coding_sequence': '', 'overhang_sequence': ''},
                '',
            ),
        ],
    )
    def test_itp_schematic(self, params, expected):
        ax = utils.itp_schematic(**params)
        # check that E P A overhang patches are drawn
        N = 3 if params.get('overhang_sequence', None) == '' else 4
        assert len(ax.patches) == N
        # checks the texts (patches, nucleotides, amino-acids)
        texts = [
            c.get_text() for c in ax.get_children() if hasattr(c, 'get_text')
        ]
        if N == 4:
            assert texts[:N] == ['protected overhang', 'E', 'P', 'A']
        else:
            assert texts[:N] == ['E', 'P', 'A']
        assert ''.join(texts[N:]) == expected

    @pytest.mark.parametrize(
        'inpt, params, expected',
        [
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
        ],
    )
    def test_dict_to_tuple(self, inpt, params, expected):
        assert utils.dict_to_tuple(inpt, **params) == expected

    @pytest.mark.parametrize(
        'inpt, params, expected',
        [
            ([], {}, []),
            (['x', 'y', 'x', 'x'], {}, ['x', 'y', 'x.1', 'x.2']),
            (['x', 'y', 'x', 'x'], {'sep': '-'}, ['x', 'y', 'x-1', 'x-2']),
            (['x', 'y', 'x.1', 'x'], {}, ['x', 'y', 'x.1', 'x.1.1']),
        ],
    )
    def test_dedup_names(self, inpt, params, expected):
        assert utils.dedup_names(inpt, **params) == expected
