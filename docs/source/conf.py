import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'tstrippy'
copyright = '2024, Salvatore Ferrone'
author = 'Salvatore Ferrone'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',"sphinx.ext.napoleon","nbsphinx"]

templates_path = ["_templates"]
exclude_patterns = ['tstrippy/tstrippy/lib/*.so','.ipynb_checkpoints/*']

master_doc = 'index'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = []

# -- Options for configuring jupyter notebooks --------------------------------
nbsphinx_prompt_width = "0"
# nbsphinx_execute = 'always' # this will execute the notebooks before building the documentation. I'm dropping it because it's too difficult to manage the dependencies of the notebooks.