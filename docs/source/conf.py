# docs/source/conf.py — cleaned for autosummary + markdown support

import sys
from pathlib import Path

# -- Path setup --------------------------------------------------------------

# Add your project src so autodoc can import it
sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

# -- Project information -----------------------------------------------------

project = "deeplotyper"
author = "Eli Niktab"
copyright = "2025, Eli Niktab"
# The full version, including alpha/beta/rc tags
release = "2025.10.2-alpha"


# -- General configuration ---------------------------------------------------

# Sphinx extensions
extensions = [
    "sphinx.ext.autodoc",           # core autodoc support
    "sphinx.ext.autosummary",       # generate stub .rst files
    "sphinx.ext.napoleon",          # Google/NumPy docstrings
    "sphinx_autodoc_typehints",     # move typehints into descriptions
    "sphinx.ext.viewcode",          # add links to source
    "sphinx.ext.todo",              # collect TODOs
    "sphinx.ext.coverage",          # coverage report
    "sphinx.ext.intersphinx",       # link to external docs
    "myst_parser",                  # Markdown support via MyST
]

# Generate autosummary pages automatically
autosummary_generate = True
autosummary_imported_members = True

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# Autodoc settings
autodoc_typehints = "description"
autodoc_member_order = "groupwise"

# Include TODOs in the output
todo_include_todos = True

# Recognize both .rst and .md
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Exclude build and system files
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Intersphinx configuration ----------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}


# -- Options for HTML output ------------------------------------------------

html_theme = "furo"
html_static_path = ["_static"]
templates_path = ["_templates"]

# Path to your logo & favicon (if you have them)
html_logo = "_static/logo.png"
html_favicon = "_static/favicon.ico"

html_theme_options = {
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
}
