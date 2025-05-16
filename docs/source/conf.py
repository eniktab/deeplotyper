# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


import sys
import subprocess, pathlib
#import sphinx_rtd_theme
sys.path.insert(
    0,
    str(pathlib.Path(__file__).parents[2] / 'src')
)


project = 'deeplotyper'
copyright = '2025, Eli Niktab'
author = 'Eli Niktab'
release = '2025.10.2-alpha'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# ---- auto-generate the API rst files on each build ----
subprocess.run([
    'sphinx-apidoc',
    '--separate',
    '--no-toc',
    '-o',
    str(pathlib.Path(__file__).parent),                      # → docs/source/
    str(pathlib.Path(__file__).parents[2] / 'src' / 'deeplotyper')  # → src/deeplotyper
], check=True)


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
    'sphinx.ext.viewcode',
    "sphinx.ext.autosummary",
    'sphinx.ext.todo',           # mark TODOs in docs
    'sphinx.ext.coverage',       # doc-coverage reports
    'sphinx.ext.intersphinx',    # link to Python stdlib, NumPy, etc.
]

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}

autosummary_generate = True
autodoc_typehints = "description"
autodoc_member_order = 'groupwise'
napoleon_google_docstring = True
napoleon_numpy_docstring = True


# Theme and Logo
# so you can get an HTML theme, e.g.:
html_theme = 'furo'
#html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_static_path = ['_static']
templates_path = ['_templates']
html_logo = "_static/logo.png"
html_favicon = "_static/favicon.ico"
html_theme_options = {
  "sidebar_hide_name": False,
  "navigation_with_keys": True,
}
