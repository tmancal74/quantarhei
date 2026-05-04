import os
import sys
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _pkg_version

sys.path.insert(0, os.path.abspath("../../"))

# -- Project information -------------------------------------------------------

project = "Quantarhei"
copyright = "2016-25, Tomas Mancal"
author = "Tomas Mancal"

try:
    version = _pkg_version("quantarhei")
except PackageNotFoundError:
    version = "unknown"
release = version

# -- General configuration -----------------------------------------------------

extensions = [
    "matplotlib.sphinxext.plot_directive",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "numpydoc",
]

templates_path = ["_templates"]
source_suffix = ".rst"
master_doc = "index"
language = "en"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
pygments_style = "sphinx"

# -- Autodoc -------------------------------------------------------------------

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
autodoc_member_order = "bysource"
autosummary_generate = True

# -- Napoleon (NumPy/Google docstrings) ----------------------------------------

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_use_admonition_for_notes = True
napoleon_use_rtype = False

# -- numpydoc ------------------------------------------------------------------

numpydoc_show_class_members = False
numpydoc_class_members_toctree = False

# -- HTML output ---------------------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_title = f"Quantarhei {release}"
html_static_path = ["_static"]
html_theme_options = {
    "navigation_with_keys": True,
    "show_nav_level": 2,
    "show_toc_level": 2,
}

# -- LaTeX output --------------------------------------------------------------

latex_documents = [
    (master_doc, "Quantarhei.tex", "Quantarhei Documentation", author, "manual"),
]

# -- Other output formats ------------------------------------------------------

man_pages = [(master_doc, "quantarhei", "Quantarhei Documentation", [author], 1)]

texinfo_documents = [
    (
        master_doc,
        "Quantarhei",
        "Quantarhei Documentation",
        author,
        "Quantarhei",
        "Molecular open quantum systems simulator.",
        "Miscellaneous",
    ),
]
