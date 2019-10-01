try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description':      'Collect seismic or infrasound waveform data and metadata from IRIS, WATC, or AVO servers or from local file sources.',
    'name':             'waveform_collection',
    'author':           'David Fee, Liam Toney, Andrew Winkelman',
    'author_email':     'dfee1@alaska.edu, ldtoney@alaska.edu, atwinkelman@alaska.edu',
    'version':          '0.1',
    'install_requires': ['obspy'],
    'scripts':          ['waveform_collection.py'],
    'data_files':       [('avo_json',['avo_json/avo_infra_calibs.json','avo_json/avo_infra_coords.json'])]
    }

setup(**config)
