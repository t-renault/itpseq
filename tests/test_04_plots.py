import os
from pathlib import Path

import pytest
from itpseq import DataSet, Replicate, Sample


class TestDataSetPlots:
    @pytest.fixture(scope='class', autouse=True)
    def setup_class(self, request, data_dir, tmp_outdir):
        request.cls.tcx_data = DataSet(
            data_dir / 'tcx_small_test', result_path=tmp_outdir
        )
        request.cls.tcx_data._clear_cache(force=True)

    def test_itp_len_plot(self):
        self.tcx_data.itp_len_plot()

    def test_itp_len_plot_row_sample(self):
        self.tcx_data.itp_len_plot(row='sample')


class TestSamplePlots:
    @pytest.fixture(scope='class', autouse=True)
    def setup_class(self, request, data_dir, tmp_outdir):
        request.cls.tcx = DataSet(
            data_dir / 'tcx_small_test', result_path=tmp_outdir
        )['nnn15.tcx']
        request.cls.tcx.dataset._clear_cache(force=True)

    def test_itp_len_plot(self):
        self.tcx.itp_len_plot()

    def test_hmap_EP(self):
        self.tcx.hmap('E', 'P')

    def test_hmap_EP_min_peptide(self):
        self.tcx.hmap('E', 'P', min_peptide=3)

    def test_hmap_pos(self):
        self.tcx.hmap_pos()

    def test_hmap_pos_EPA(self):
        self.tcx.hmap_pos(('E', 'P', 'A'))

    def test_hmap_grid(self):
        self.tcx.hmap_grid()

    def test_volcano(self):
        self.tcx.volcano('E:P')

    def test_itoeprint(self):
        self.tcx.itoeprint(show_range=True)

    def test_itoeprint_shades(self):
        self.tcx.itoeprint('shades', show_range=True)

    def test_subset_logo(self):
        for logo_type in ['raw_freq', 'extra_counts', 'sum_log2FC']:
            self.tcx.subset_logo(
                'E:P', logo_type=logo_type, query='log2FoldChange > 0'
            )
            self.tcx.subset_logo(
                'E:P',
                logo_type=f'{logo_type}_bits',
                query='log2FoldChange > 0',
            )

    def test_logo(self):
        self.tcx.logo()

    def test_all_logos(self):
        self.tcx.all_logos()
