import filecmp
from itertools import product
from pathlib import Path
from unittest.mock import Mock

import itpseq.__main__ as cli
import pytest
from click.testing import CliRunner


class TestCLI:
    @pytest.mark.parametrize(
        'inpt, expected',
        [
            ('1', (1, 1)),
            ('1-', (1, 0)),
            ('-10', (0, 10)),
            ('1-10', (1, 10)),
        ],
    )
    def test_cli_arg_range(self, inpt, expected):
        assert cli.arg_range(None, None, inpt) == expected

    def test_cli_parse(self, data_dir, tmp_outdir):
        path = data_dir / 'tcx_small_test'
        test_files = list(path.glob('*.fastq'))

        runner = CliRunner()
        result = runner.invoke(
            cli.parse, [*map(str, test_files), '-o', str(tmp_outdir)]
        )

        assert (
            result.exit_code == 0
        ), f'"parse" failed with error: {result.output}'

        files = [
            ''.join(p)
            for p in product(
                (f.name.removesuffix('.assembled.fastq') for f in test_files),
                ('.aa.itp.txt', '.nuc.itp.txt', '.itp.log', '.itp.json'),
            )
        ]

        ok, nok, err = filecmp.cmpfiles(path, tmp_outdir, files, shallow=False)
        assert nok == [] and err == []

    def test_cli_parse_a1(self, data_dir, tmp_outdir):
        path = data_dir / 'tcx_small_test'
        res_path = data_dir / 'fake_output_data' / 'tcx_a1_ATG'
        test_file = 'nnn15_noa1.assembled.fastq'
        runner = CliRunner()
        result = runner.invoke(
            cli.parse,
            [
                '-a1',
                'GTATAAGGAGGAAAAAATATG',
                str(path / test_file),
                '-o',
                str(tmp_outdir),
            ],
        )
        assert (
            result.exit_code == 0
        ), f'"parse" failed with error: {result.output}'

        files = [
            ''.join(p)
            for p in product(
                [test_file.removesuffix('.assembled.fastq')],
                ('.aa.itp.txt', '.nuc.itp.txt', '.itp.log', '.itp.json'),
            )
        ]

        ok, nok, err = filecmp.cmpfiles(
            res_path, tmp_outdir, files, shallow=False
        )
        assert nok == [] and err == []

    def test_cli_format(self, data_dir):
        path = data_dir / 'tcx_small_test'
        runner = CliRunner()
        result = runner.invoke(
            cli.format, ['--limit', '3', str(path / 'nnn15_noa1.nuc.itp.txt')]
        )
        assert (
            result.exit_code == 0
        ), f'"format" failed with error: {result.output}'
        assert result.stdout.splitlines() == [
            '#                    [E][P][A]               ',
            '                     ATGGGACGC cccgcagtatct  ',
            '      ATGAGTTACAAAGGCAACTCGGAA caggtagcatatc ',
            '                     ATGGAAGAG gcccatgccattcc',
        ]

    def test_cli_report(self, data_dir, tmp_outdir):
        path = data_dir / 'tcx_small_test'
        runner = CliRunner()
        result = runner.invoke(cli.report, [str(path), '-o', str(tmp_outdir)])
        assert (
            result.exit_code == 0
        ), f'"parse" failed with error: {result.output}'

    def test_cli_help(self, monkeypatch):
        runner = CliRunner()
        result = runner.invoke(cli.help, ['--help'])
        assert (
            result.exit_code == 0
        ), f'"help --help" failed with error: {result.output}'
        mock_open = Mock()
        monkeypatch.setattr('webbrowser.open', mock_open)
        result = runner.invoke(cli.help)
        assert (
            result.exit_code == 0
        ), f'"help" failed with error: {result.output}'
        mock_open.assert_called_once()

    def test_cli_completion(self):
        runner = CliRunner()
        result = runner.invoke(cli.completion, ['bash'])
        assert (
            result.exit_code == 0
        ), f'"help --help" failed with error: {result.output}'
