.. _command_line:

=====================
Command line analysis
=====================

After :ref:`installation<installation>`, the ``itpseq`` command should be
available in the shell.

.. _cli_parsing:

Parsing the fastq files
-----------------------

A parsing step is required to extract the inverse toeprints from the fastq
files. This operation will create four file for each input fastq file
``<file_prefix>.fastq`` (or ``<file_prefix>.assembled.fastq``):

* the inverse-toeprint sequences as nucleotides (``<file_prefix>.nuc.txt``) 

.. code-block:: console

  #                    [E][P][A]              
                       ATGGGACGCCCCGCAGTATCT  
        ATGAGTTACAAAGGCAACTCGGAACAGGTAGCATATC 
                       ATGGAAGAGGCCCATGCCATTCC
                          ATGAATCGAAACATGTTT  
  ATGACTATGTTTCTTGGACACACATAAGGGAACTAGTTAGGG  
                       ATGCTATAATAGGTCAAGCACCA
                 ATGACCAATCCGTAGGACTAACGCCACAT
           ATGTAGCCGGGCAAGGAGATCCGCACCTCGCGC  
                          ATGTAACTATACGACGTCG 

* the inverse-toeprint sequences as amino-acids (``<file_prefix>.aa.itp.txt``)

.. code-block:: console

  #      EPA
         mGR
    mSYKGNSE
         mEE
          mN
  mTMFLGHT*G
         mL*
       mTNP*
     m*PGKEI
          m*

* metadata as JSON (``<file_prefix>.itp.json``)
* a log file (``<file_prefix>.itp.log``)

The ``itpseq parse`` command processes the specified FASTQ files and generates  
the four corresponding output files for each input file:

.. code-block:: console

    itpseq parse *.assembled.fastq

Options are available to specify the output directory, parameters to filter the
reads, or the adaptors that were used in the design. For more details run:

.. code-block:: console

    itpseq parse --help

.. _cli_report:

Generating the report
---------------------

A report with a default set of analyses and graphs can be obtained using:

.. code-block:: console

    itpseq report data_dir

Where ``data_dir`` is the directory containing the data files. For the current
directory, run: 

.. code-block:: console

    itpseq report .

