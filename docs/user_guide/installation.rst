.. _installation:

============
Installation
============

The code for `itpseq` is available on GitHub
(`<https://github.com/t-renault/itpseq>`__) under an GPL-3.0 open source license.

We assume that `Python <https://www.python.org/downloads/>`_ and the `pip
<https://packaging.python.org/en/latest/tutorials/installing-packages/>`_
package manager are already installed on your machine.

The `itpseq` module and its command line tools can be installed/updated from `PyPI <https://pypi.org/project/itpseq/>`_ by running the following command:

.. code-block:: console

    pip install -U itpseq

This should make the `itpseq` command available in your shell, as well as a python module.

In order to check the correct installation, run:

.. code-block:: console

    itpseq --help

which should output the help page:

.. code-block:: console

    Usage: itpseq [OPTIONS] COMMAND [ARGS]...

      itpseq - A command-line tool for the analysis of iTP-Seq data.

    Options:
      -v, --version  Show the version and exit.
      -h, --help     Show this message and exit.

    Commands:
      help    Open the HTML documentation in the default web browser.
      parse   Parse and filter the assembled iTP-Seq FASTQ files to produce...
      report  Generate a report for the dataset in the specified DIRECTORY.

      To get help on a specific command, run: itpseq COMMAND --help

From within a python interpreter or python script:

.. ipython:: python

    import itpseq

    print(itpseq.__version__)
