# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys

sys.path.insert(0, os.path.abspath("../../"))
from sectionproperties import __version__ as ver

# -- Project information -----------------------------------------------------

project = "sectionproperties"
copyright = "2022, Robbie van Leeuwen"
author = "Robbie van Leeuwen"

# The short Major.Minor.Build version
_v = ver.split(".")
_build = "".join([c for c in _v[2] if c.isdigit()])
version = _v[0] + "." + _v[1] + "." + _build
# The full version, including alpha/beta/rc tags
release = ver


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
    "matplotlib.sphinxext.plot_directive",
]

autodoc_member_order = "bysource"
autosummary_generate = True  # Turn on sphinx.ext.autosummary
autosummary_mock_imports = []
autoclass_content = "both"  # Add __init__ doc (ie. params) to class summaries
html_show_sourcelink = (
    False  # Remove 'view source code' from top of page (for html, not python)
)
autodoc_inherit_docstrings = True  # If no docstring, inherit from base class
add_module_names = False  # Remove namespaces from class/method signatures
# nbsphinx_allow_errors = True  # Continue through Jupyter errors
# nbsphinx_execute_arguments = [
#     "--InlineBackend.figure_formats={'svg', 'pdf'}",
#     "--InlineBackend.rc=figure.dpi=96",
# ]
typehints_use_rtype = False  # document return type as part of the :return: directive
typehints_defaults = "comma"  # adds a default annotation after the type

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# intersphinx mapping
intersphinx_mapping = {
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"

# # Define the json_url for our version switcher.
# json_url = (
#     "https://robbievanleeuwen.github.io/concrete-properties/_static/switcher.json"
# )

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/robbievanleeuwen/section-properties",
            "icon": "fab fa-github-square",
            "type": "fontawesome",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/sectionproperties/",
            "icon": "fas fa-box",
            "type": "fontawesome",
        },
    ],
    "use_edit_page_button": True,
    "logo": {
        "image_light": "logo.png",
        "image_dark": "logo_dark.png",
    },
    "navbar_start": ["navbar-logo"],
    # "navbar_start": ["navbar-logo", "version-switcher"],
    # "switcher": {
    #     "json_url": json_url,
    #     # "version_match": "latest",
    #     "version_match": release,
    # },
}

html_context = {
    "github_user": "robbievanleeuwen",
    "github_repo": "section-properties",
    "github_version": "master",
    "doc_path": "docs/source/",
    "default_mode": "light",
}