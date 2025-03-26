# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import os
import sys

sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from itpseq import __version__

project = 'iTP-Seq'
copyright = '2025, Thibaud Renault'
author = 'Thibaud Renault'
release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'numpydoc',
    'sphinx_design',
    'nbsphinx',
    'sphinx_gallery.load_style',
    'sphinx.ext.linkcode',
    'sphinx.ext.doctest',
    'IPython.sphinxext.ipython_console_highlighting',
    'IPython.sphinxext.ipython_directive',
    #'myst_nb',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
source_suffix = ['.rst', '.md']
nbsphinx_execute = 'never'

html_show_sourcelink = False

autosummary_generate = True
autodoc_member_order = 'groupwise'
add_function_parentheses = False

nbsphinx_thumbnails = {
    'notebooks/itoeprint_tutorial': '_static/itoeprint_tutorial.png',
    #'notebooks/report_tutorial': '_static/report_tutorial.png',
    'notebooks/report_tutorial': '_static/report_page.png',
}

import inspect
import warnings

import itpseq


def linkcode_resolve(domain, info) -> str | None:
    """
    Determine the URL corresponding to Python object
    """
    if domain != 'py':
        return None

    modname = info['module']
    fullname = info['fullname']

    submod = sys.modules.get(modname)
    if submod is None:
        return None

    obj = submod
    for part in fullname.split('.'):
        try:
            with warnings.catch_warnings():
                # Accessing deprecated objects will generate noisy warnings
                warnings.simplefilter('ignore', FutureWarning)
                obj = getattr(obj, part)
        except AttributeError:
            return None

    try:
        fn = inspect.getsourcefile(inspect.unwrap(obj))
    except TypeError:
        try:  # property
            fn = inspect.getsourcefile(inspect.unwrap(obj.fget))
        except (AttributeError, TypeError):
            fn = None
    if not fn:
        return None

    try:
        source, lineno = inspect.getsourcelines(obj)
    except TypeError:
        try:  # property
            source, lineno = inspect.getsourcelines(obj.fget)
        except (AttributeError, TypeError):
            lineno = None
    except OSError:
        lineno = None

    if lineno:
        linespec = f'#L{lineno}-L{lineno + len(source) - 1}'
    else:
        linespec = ''

    fn = os.path.relpath(fn, start=os.path.dirname(itpseq.__file__))

    if '+' in itpseq.__version__:
        return f'https://github.com/t-renault/itpseq/blob/main/itpseq/{fn}{linespec}'
    else:
        return (
            f'https://github.com/t-renault/itpseq/blob/'
            f'{itpseq.__version__}/itpseq/{fn}{linespec}'
        )


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'classic'
html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_css_files = ['css/itpseq.css']
