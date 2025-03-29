itpseq: A Python module and command-line tool for the analysis of iTP-Seq data.
==================================

[![PyPI Version](https://img.shields.io/pypi/v/itpseq.svg)](https://pypi.org/project/itpseq/)
[![License](https://img.shields.io/pypi/l/itpseq.svg)](https://github.com/t-renault/itpseq/blob/main/LICENSE)
[![Tests](https://github.com/t-renault/itpseq/actions/workflows/cross-platform-testing.yml/badge.svg)](https://github.com/t-renault/itpseq/actions/workflows/cross-platform-testing.yml)
[![Code Coverage](https://codecov.io/gh/t-renault/itpseq/branch/main/graph/badge.svg)](https://codecov.io/gh/t-renault/itpseq)
[![Downloads](https://pepy.tech/badge/itpseq)](https://pepy.tech/project/itpseq)

The latest documentation is available at [itpseq.readthedocs.io](https://itpseq.readthedocs.io).

Dependencies
------------
`itpseq` runs on Python 3.9+ on Linux, MacOS, and Windows.

Installation requires [regex](https://pypi.org/project/regex/),
[numpy](https://numpy.org/),
[pandas](https://pandas.pydata.org/),
[matplotlib](https://matplotlib.org/),
[seaborn](https://seaborn.pydata.org/),
[pydeseq2](https://pydeseq2.readthedocs.io/en/stable/),
[logomaker](https://logomaker.readthedocs.io/),
[Jinja2](https://jinja.palletsprojects.com/en/stable/),
[weasyprint](https://pypi.org/project/weasyprint/),
[ipython](https://ipython.org/),
[click](https://click.palletsprojects.com/),
and [ipywidgets](https://ipywidgets.readthedocs.io/en/stable/).
The dependencies are installed automatically when installing from pip.


Installation
------------
The latest stable release (and required dependencies) can be installed from PyPI:

    pip install itpseq

To update the package to the latest version, run:

    pip install -U itpseq

After installation, if Python and pip are configured correctly, the "itpseq" command should be available in the terminal. If not, one should ensure the location of the executable is in the `PATH`. 

Citing
------


Known issues
------------
On MacOS and Windows, the "weasyprint" library 	that is used to produce the PDF version of the report seems to have an issue and exporting a PDF report is currently not possible.
It is nevertheless possible to produce an HTML report by running:

    itpseq report -o report.html data_directory 
  
in place of:

    itpseq report data_directory 
