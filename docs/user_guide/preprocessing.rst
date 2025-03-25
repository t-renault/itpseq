.. _preprocessing:

========================================
Pre-processing data to use with `itpseq`
========================================

Pre-requisite
-------------

The pre-requisite to use ``itpseq`` is to have **one** fastq file per replicate.

The `fastq <https://en.wikipedia.org/wiki/FASTQ_format>`_ files are text files
with four lines per read:

.. code-block:: console

    @E00490:646:H7HGCCCX2:8:1101:5203:2206 1:N:0:CACGAT
    GTATAAGGAGGAAAAAATATGCCAGCATCCGTCTAGTCCACGCGGTATCTCGGTGTGACTGACTGAGAATTTCTGTAGGCACCATCAAT
    +
    \V\MaMDRe``jVeVo``JJV\`VRNVNM9\eHf`RIRjkaV\\fjoo\eoaojjssoVosVVkaaVsssssVsfojossee`j\jeeM

If paired-end sequencing was performed, it is necessary to :ref:`merge the files<pairing>`.
To be able to :ref:`automatically import<auto_loading>` the files in a
:class:`DataSet`, it is best to follow a :ref:`naming convention<naming>`.


.. _pairing:

Merging paired-end sequencing files
-----------------------------------

In case of paired-end sequencing, it is first necessary to pair the two reads.
This can be done using tools such as `PEAR
<https://cme.h-its.org/exelixis/web/software/pear/>`_ or `Pandaseq
<https://pandaseq-tutorial.readthedocs.io/en/latest/tutorial/>`_.

For instance, with PEAR, assuming a forward ``nnn15_noa1_1.fastq`` and reverse
``nnn15_noa1_2.fastq``, those can be merged into a single
``nnn15_noa1.assembled.fastq`` using:

.. code-block:: console

    pear -u0 -c126 -f nnn15_noa1_1.fastq -r nnn15_noa1_2.fastq -o nnn15_noa1

After pairing of the reads, a dataset with a control (``noa``) and an
antibiotic-treated condition (``tcx``), with 3 replicates each might look like:

.. code-block:: console

    nnn15_noa1.assembled.fastq
    nnn15_noa2.assembled.fastq
    nnn15_noa3.assembled.fastq
    nnn15_tcx1.assembled.fastq
    nnn15_tcx2.assembled.fastq
    nnn15_tcx3.assembled.fastq

.. _naming:

Naming conventions
------------------

``itpseq`` relies on filenames to :ref:`automatically load<auto_loading>` data files and match the
references. Although not strictly required we recommend to use filenames without
spaces and with short descriptors of the sample properties and groups.

For instance, we use ``nnn15`` as a prefix to indicate that a random NNN15
library was used. In case of an antibiotic treatment vs no-antibiotic control we
use ``noa`` for the control and a short identifier such as ``tcx`` for the
antibiotic-treated samples. We use a number to indicate the replicate. A
filename for the second replicate of an erythromycin-treated sample, with a
NNN15 library might be named ``nnn15_ery2.assembled.fastq``.

During the analysis, ``itpseq`` can combine automatically the replicates into a
``Sample``, and assign a reference to it. See :ref:`loading_data` for more
details.
