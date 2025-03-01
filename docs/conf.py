# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import sys, os
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'iTP-Seq'
copyright = '2025, Thibaud Renault'
author = 'Thibaud Renault'
release = '0.0.1a3'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'numpydoc',
    'nbsphinx',
    'sphinx_gallery.load_style',
    #'myst_nb',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
source_suffix = ['.rst', '.md']
nbsphinx_execute = 'never'

autosummary_generate = True
autodoc_member_order = 'groupwise'

nbsphinx_thumbnails = {
    'notebooks/itoeprint_tutorial': '_static/itoeprint_tutorial.png',
    #'notebooks/report_tutorial': '_static/report_tutorial.png',
    'notebooks/report_tutorial': '_static/report_page.png',
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
