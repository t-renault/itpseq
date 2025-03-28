import os
from pathlib import Path

import pytest

from itpseq import DataSet


class TestReport:
    @pytest.fixture(scope='class', autouse=True)
    def setup_class(self, request, data_dir, tmp_outdir):
        request.cls.tcx_data = DataSet(
            data_dir / 'tcx_small_test', result_path=tmp_outdir
        )

    def test_report_html(self, tmp_outdir):
        self.tcx_data._clear_cache(force=True)
        self.tcx_data.report(output=Path(tmp_outdir) / 'report.html')
