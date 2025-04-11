.. _user_guide:

==========
User Guide
==========


Welcome to the iTP-Seq user guide. This documentation will help you understand the iTP-Seq method and learn how to analyze iTP-Seq data using the ``itpseq`` library.

.. _start:

Getting started
---------------

iTP-Seq (inverse toeprinting coupled to next-generation sequencing) is a powerful method for characterizing bacterial translation landscapes.
This guide covers the method and data analysis workflow.

* :ref:`whatis_itpseq` - Overview of the iTP-Seq method, applications, and comparison to other techniques
* :ref:`installation` - Step-by-step instructions for installing the ``itpseq`` software package
* :ref:`preprocessing` - Guidelines for organizing, formatting, and preparing data for analysis
* :ref:`command_line` - Documentation for the ``itpseq`` command line interface
* :ref:`loading_data` - Instructions for loading and manipulating iTP-Seq data in Python
* :ref:`tutorials` - Practical examples and walkthroughs for common analysis tasks

Analysis approaches
-------------------

The ``itpseq`` package can be used in two ways:

1. **Command line interface** - Suitable for standardized analyses and batch processing.
See the :ref:`command_line` section for details.

2. **Python module** - For custom analysis pipelines and interactive exploration.
Refer to the :ref:`loading_data` section for getting started and :ref:`api` for the Python API documentation.

Need help?
----------

For new users, we recommend beginning with   :ref:`whatis_itpseq` to understand the method, followed by the :ref:`installation` guide to set up the software.
The :ref:`tutorials` section demonstrates typical workflows.

For bug reports or feature requests, please use the GitHub issue tracker: https://github.com/t-renault/itpseq/issues


.. toctree::
    :maxdepth: 2
    :hidden:
    
    whatis_itpseq
    installation
    preprocessing
    command_line
    loading_data
    tutorials
    
