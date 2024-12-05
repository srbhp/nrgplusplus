import os
import sys

sys.path.insert(0, os.path.abspath("."))
# -- Project Information----------------------------------
project = "nrgplusplus"
copyright = "2017-2023, Saurabh Pradhan"
author = "Saurabh Pradhan"


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = "4.5"


# The `extensions` list should already be in here from `sphinx-quickstart`
extensions = [
    # there may be others here already, e.g. 'sphinx.ext.mathjax'
    "breathe",
    "exhale",
]

# Setup the breathe extension
breathe_projects = {"nrgplusplus": "../build/_doxygen/xml"}
breathe_default_project = "nrgplusplus"

# Setup the exhale extension
exhale_args = {
    # These arguments are required
    "containmentFolder": "./api",
    "rootFileName": "library_root.rst",
    "doxygenStripFromPath": "../nrgcore/include",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle": "Library API",
    # Suggested optional arguments
    "createTreeView": True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin": "INPUT = ../nrgcore/include",
    "pageLevelConfigMeta": ":github_url: https://github.com/srbhp/nrgplusplus",
}

# Tell sphinx what the primary language being documented is.
primary_domain = "cpp"

# Tell sphinx what the pygments highlight language should be.
highlight_language = "cpp"

# The short X.Y version.
# import exhale
# NOTE: this is the companion site for Exhale, which is why I'm setting the
#       version to be the same.  For your own projects, you would NOT do this!
version = "0.0.1"  # exhale.__version__
# The full version, including alpha/beta/rc tags.
# release = "rc"  # exhale.__version__

# -- Options for HTML output -------------------------------------------------

# [[[ begin theme marker ]]]
# The name of the Pygments (syntax highlighting) style to use.
# `sphinx` works very well with the RTD theme, but you can always change it
pygments_style = "sphinx"

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"
# sidebar
# html_sidebars = {
#    "**": [
#        "about.html",
#        "navigation.html",
#        "relations.html",
#        "searchbox.html",
#        "donate.html",
#    ]
# }
html_theme_options = {
    "header_links_before_dropdown": 4,
    "icon_links": [
        {
            # Label for this link
            "name": "GitHub",
            # URL where the link will redirect
            "url": "https://github.com/srbhp/nrgplusplus",  # required
            # Icon class (if "type": "fontawesome"), or path to local image (if "type": "local")
            "icon": "fa-brands fa-square-github",
            # The type of image to be used (see below for details)
            "type": "fontawesome",
        }
    ],
}
html_theme_options = {}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
# [[[ end theme marker ]]]

rst_epilog = ".. |theme| replace:: ``{0}``".format(html_theme)
