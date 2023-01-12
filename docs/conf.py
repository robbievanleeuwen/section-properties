"""Sphinx configuration."""

# project information
project = "sectionproperties"
author = "Robbie van Leeuwen"
copyright = "2023, Robbie van Leeuwen"

# sphinx config
add_module_names = False  # remove namespaces from class/method signatures
templates_path = ["_templates"]

# extensions
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    # "IPython.sphinxext.ipython_console_highlighting",
    # "matplotlib.sphinxext.plot_directive",
    "myst_parser",
    # "nbsphinx",
    # "sphinx_autodoc_typehints",
    "sphinx_click",
    "sphinx_copybutton",
    "sphinxext.opengraph",
]

# extension config
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autoclass_content = "both"  # add __init__ doc (ie. params) to class summaries
autodoc_inherit_docstrings = True  # If no docstring, inherit from base class
add_module_names = False  # Remove namespaces from class/method signatures
# nbsphinx_execute_arguments = [
#     "--InlineBackend.figure_formats={'svg', 'pdf'}",
#     "--InlineBackend.rc=figure.dpi=96",
# ]
# typehints_use_rtype = False  # document return type as part of the :return: directive
# typehints_defaults = "comma"  # adds a default annotation after the type

# intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "shapely": ("https://shapely.readthedocs.io/en/stable/", None),
}

# html theme
html_theme = "furo"
html_static_path = ["_static"]
html_theme_options = {
    "light_logo": "logo-light-mode.png",  # add light mode logo
    "dark_logo": "logo-dark-mode.png",  # add dark mode logo
    "sidebar_hide_name": True,  # hide name of project in sidebar (already in logo)
    "source_repository": "https://github.com/robbievanleeuwen/section-properties",
    "source_branch": "master",
    "source_directory": "docs/",
}
pygments_style = "sphinx"
pygments_dark_style = "monokai"
