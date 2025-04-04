import filecmp
from itertools import product
from pathlib import Path

import itpseq.__main__ as cli
import pytest
from click.testing import CliRunner


class TestCLI:
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
        if nok:
            print(f'Files are not equal: {nok}')
        if err:
            print(f'Files could not be read: {err}')
        assert nok == [] and err == []

    def test_cli_report(self, data_dir, tmp_outdir):
        path = data_dir / 'tcx_small_test'
        runner = CliRunner()
        result = runner.invoke(cli.report, [str(path), '-o', str(tmp_outdir)])
        assert (
            result.exit_code == 0
        ), f'"parse" failed with error: {result.output}'

    def test_cli_help(self):
        runner = CliRunner()
        result = runner.invoke(cli.help, ['--help'])
        assert (
            result.exit_code == 0
        ), f'"help --help" failed with error: {result.output}'
        result = runner.invoke(cli.help)
        assert (
            result.exit_code == 0
        ), f'"help" failed with error: {result.output}'

    def test_cli_completion(self):
        runner = CliRunner()
        result = runner.invoke(cli.completion, ['bash'])
        assert (
            result.exit_code == 0
        ), f'"help --help" failed with error: {result.output}'
