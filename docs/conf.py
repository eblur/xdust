# Configuration file for the Sphinx documentation builder.
import os
import sys

# Make the xdust package importable from the src layout
sys.path.insert(0, os.path.abspath('../src'))

# -- Project information -----------------------------------------------------
project = 'xdust'
copyright = '2024, Lia Corrales'
author = 'Lia Corrales'
version = '0.1'
release = '0.1'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    # 'sphinx_autodoc_typehints', # requires a different format: https://pypi.org/project/sphinx-autodoc-typehints/
]

autosummary_generate = True

autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
}

templates_path = ['_templates']
exclude_patterns = ['_build']

master_doc = 'index'
source_suffix = '.rst'
language = 'en'

# -- Options for HTML output -------------------------------------------------
html_theme = 'furo'

html_theme_options = {}

html_static_path = []
