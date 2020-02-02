import os
import sys
sys.path.insert(0, os.path.abspath('../waveform_collection'))

project = 'waveform_collection'

html_show_copyright = False

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.intersphinx',
              'recommonmark',
              'sphinx.ext.viewcode',
              'sphinxcontrib.apidoc',
              'sphinx.ext.mathjax'
              ]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'

napoleon_numpy_docstring = False

master_doc = 'index'

autodoc_mock_imports = ['numpy',
                        'obspy'
                        ]

apidoc_module_dir = '../waveform_collection'

apidoc_output_dir = 'api'

apidoc_separate_modules = True

apidoc_toc_file = False

# Exclude WATC-related stuff from docs for now
apidoc_excluded_paths = ['local/cd11.py',
                         'local/clf.py',
                         'local/css.py',
                         'local/smart24.py']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://docs.scipy.org/doc/numpy', None),
    'obspy': ('https://docs.obspy.org/', None)
}
