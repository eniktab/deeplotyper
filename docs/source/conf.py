# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
#import sphinx_rtd_theme
sys.path.insert(0, os.path.abspath('../../src'))


project = 'deeplotyper'
copyright = '2025, Eli Niktab'
author = 'Eli Niktab'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
    'sphinx.ext.viewcode',
]

autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# so you can get an HTML theme, e.g.:
html_theme = 'alabaster'
#html_theme = 'sphinx_rtd_theme'
#html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
