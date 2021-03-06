# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import dask.dataframe
import xarray

sys.path.insert(0, os.path.abspath('../src/local_pcangsd'))


# -- Project information -----------------------------------------------------

project = "local_pcangsd"
copyright = "2022, Alexis Simon"
author = "Alexis Simon"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx.ext.autosummary",
    "sphinx.ext.githubpages",
    "sphinx.ext.intersphinx",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# FIXME: Workaround for linking xarray module
# For some reason, intersphinx is not able to to link xarray objects.
# https://github.com/pydata/xarray/issues/4279
xarray.Dataset.__module__ = "xarray"
xarray.DataArray.__module__ = "xarray"

dask.dataframe.DataFrame.__module__ = "dask.dataframe"
intersphinx_mapping = dict(
    dask=("https://docs.dask.org/en/stable/", None),
    xarray=("https://xarray.pydata.org/en/stable/", None),
    zarr=("https://zarr.readthedocs.io/en/stable", None),
    numpy=("https://numpy.org/doc/stable/", None),
    python=("https://docs.python.org/3", None),
    # add sgkit
)

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"


html_theme_options = {
    "repository_url": "https://github.com/alxsimon/local_pcangsd",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_download_button": False,
    "use_fullscreen_button": False,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
