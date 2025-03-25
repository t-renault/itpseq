.. _loading_data:

==============================
Loading iTP-Seq data in Python
==============================

.. _auto_loading:

Automatic loading from a directory
----------------------------------

The easiest approach to create a :class:`DataSet` is to use a consistent format of
the file names (see :ref:`naming`).

The :ref:`parsing step<cli_parsing>` creates four file for each input fast files:

* the inverse-toeprint sequences as nucleotides (``<file_prefix>.nuc.txt``) 
* the inverse-toeprint sequences as amino-acids (``<file_prefix>.aa.itp.txt``)
* metadata as JSON (``<file_prefix>.itp.json``)
* a log file (``<file_prefix>.itp.log``)

All the files share the same prefix and the JSON files are used to identify the
replicates.

Default behavior
~~~~~~~~~~~~~~~~

By default :class:`DataSet` expects a prefix with the ``XXX_YYYDD`` format. ``XXX``
(alphanumeric) will be assigned as a "lib-type" key, ``YYY`` (letters) as a "sample"
key, and ``DD`` (digits) as the "replicate". For example ``nnn15_noa1``.

Therefore a directory containing 3 "noa" and 3 "tcx" replicates would look like:

.. code-block:: console

    nnn15_noa1.aa.itp.txt       nnn15_noa3.aa.itp.txt       nnn15_tcx2.aa.itp.txt
    nnn15_noa1.assembled.fastq  nnn15_noa3.assembled.fastq  nnn15_tcx2.assembled.fastq
    nnn15_noa1.itp.json         nnn15_noa3.itp.json         nnn15_tcx2.itp.json
    nnn15_noa1.itp.log          nnn15_noa3.itp.log          nnn15_tcx2.itp.log
    nnn15_noa1.nuc.itp.txt      nnn15_noa3.nuc.itp.txt      nnn15_tcx2.nuc.itp.txt
    nnn15_noa2.aa.itp.txt       nnn15_tcx1.aa.itp.txt       nnn15_tcx3.aa.itp.txt
    nnn15_noa2.assembled.fastq  nnn15_tcx1.assembled.fastq  nnn15_tcx3.assembled.fastq
    nnn15_noa2.itp.json         nnn15_tcx1.itp.json         nnn15_tcx3.itp.json
    nnn15_noa2.itp.log          nnn15_tcx1.itp.log          nnn15_tcx3.itp.log
    nnn15_noa2.nuc.itp.txt      nnn15_tcx1.nuc.itp.txt      nnn15_tcx3.nuc.itp.txt

Loading this directory will automatically assign the 3 Replicates to 2 Samples
("tcx" and "noa"). In addition, if a sample is named "noa", it is automatically
assigned as a reference to the other samples that share the same keys (other
than "sample"):

.. ipython:: python
   :suppress:

   # save documentation directory
   # change to test data directory
   import os
   from pathlib import Path
   CWD = Path(os.getcwd())
   TEST_DIR = CWD.parent / 'tests/test_data/'
   os.chdir(TEST_DIR / 'tcx_small_test/')

.. ipython:: python

  from itpseq import DataSet
  data = DataSet('.') # current directory containing the data files
  data
  data.samples
  data.replicates

This detection is due to the default regular expression ``file_pattern``:
``(?P<lib_type>[^_]+)_(?P<sample>[^_\\d]+)(?P<replicate>\\d+)``.

The ``lib_type`` and ``sample`` keys are automatically used to group the
Replicates into a Sample and to create the Sample name.

It is possible to specify the keys to use to group the Replicates:

.. ipython:: python
   :suppress:

   os.chdir(TEST_DIR / 'tcx_small_test/')

.. ipython:: python

  DataSet('.', keys=['sample'])  # ignoring "lib_type"


Custom prefix and keys
~~~~~~~~~~~~~~~~~~~~~~

Let's imagine a dataset with two drugs (drugA and drugB), one control (noa) and
a few different concentrations for the drugs (10, 20, 30µM):

.. code-block:: console

   drugA1_10µM.itp.json  drugA3_20µM.itp.json  drugB2_30µM.itp.json
   drugA1_20µM.itp.json  drugA3_30µM.itp.json  drugB3_10µM.itp.json
   drugA1_30µM.itp.json  drugB1_10µM.itp.json  drugB3_20µM.itp.json
   drugA2_10µM.itp.json  drugB1_20µM.itp.json  drugB3_30µM.itp.json
   drugA2_20µM.itp.json  drugB1_30µM.itp.json  noa1.itp.json
   drugA2_30µM.itp.json  drugB2_10µM.itp.json  noa2.itp.json
   drugA3_10µM.itp.json  drugB2_20µM.itp.json  noa3.itp.json

The different parts of the filename can be defined through ``file_pattern``:

``(?P<sample>[^_]+)(?P<replicate>\\d+)(_(?P<concentration>\\d+µM))?``

* ``(?P<sample>[^_]+)``: match the sample name (anything but ``_``)
* ``(?P<replicate>\d+)``: match digits defining the replicate number
* ``(_(?P<concentration>\d+µM))?``: optionally match ``_`` followed by a concentration

.. ipython:: python
   :suppress:

   os.chdir(TEST_DIR / 'loading_concentrations/')

.. ipython:: python

  from itpseq import DataSet
  data = DataSet('.', file_pattern=r'(?P<sample>[^_]+)(?P<replicate>\d+)(_(?P<concentration>\d+µM))?')
  data
  data.samples

It is also possible to define the keys that will be used to assign the
replicate. For instance, using ``ref_labels={'sample': 'drugA'}`` would define
``drugA`` as a reference to the samples that match the other same keys.

.. ipython:: python

  data = DataSet('.',
                 file_pattern=r'(?P<sample>[^_]+)(?P<replicate>\d+)(_(?P<concentration>\d+µM))?',
                 ref_labels={'sample': 'drugA'},
                 )
  data


.. _manual_loading:

Manual loading
--------------

From ``Sample``/``Replicate`` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is also possible to create :class:`Replicate`, :class:`Sample`, and
:class:`DataSet` objects manually. 

.. ipython:: python
   :suppress:

   os.chdir(TEST_DIR / 'tcx_small_test/')

.. ipython:: python

  from itpseq import DataSet, Sample, Replicate

  R1 = Replicate(replicate='1', file_prefix='nnn15_tcx1') # relative to current directory
  R2 = Replicate(replicate='2', file_prefix='nnn15_tcx2')
  R3 = Replicate(replicate='3', file_prefix='nnn15_tcx3')
  N1 = Replicate(replicate='1', file_prefix='nnn15_noa1')
  N2 = Replicate(replicate='2', file_prefix='nnn15_noa2')
  N3 = Replicate(replicate='3', file_prefix='nnn15_noa3')
  
  S = Sample(replicates=[R1, R2, R3],
             name='tcx',
             reference=Sample(replicates=[N1, N2, N3], name='noa'),
            )
  S

From a dictionary
~~~~~~~~~~~~~~~~~

From a dictionary of samples/replicates:

.. ipython:: python

  data = DataSet({'tcx': [{'file_prefix': 'nnn15_tcx1'},
                          {'file_prefix': 'nnn15_tcx2'},
                          {'file_prefix': 'nnn15_tcx3'}
                         ],
                  'noa': [{'file_prefix': 'nnn15_noa1'},
                          {'file_prefix': 'nnn15_noa2', 'replicate': 'custom_name'},
                          {'file_prefix': 'nnn15_noa3'}
                         ]},
                 ref_mapping={'tcx': 'noa'})

  data

From config files (e.g. JSON)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From a JSON file (e.g. ``samples.json``):

.. code-block:: console

    {
     "tcx":[
      {"file_prefix":"nnn15_tcx1"},
      {"file_prefix":"nnn15_tcx2"},
      {"file_prefix":"nnn15_tcx3"}
     ],
     "noa":[
      {"file_prefix":"nnn15_noa1"},
      {"file_prefix":"nnn15_noa2","replicate":"custom_name"},
      {"file_prefix":"nnn15_noa3"}
     ]
    }

.. ipython:: python

  import json
  from itpseq import DataSet

  with open('samples.json') as f:
      data = DataSet(json.load(f))

  data

.. ipython:: python
   :suppress:

   # restore documentation directory
   os.chdir(CWD)
