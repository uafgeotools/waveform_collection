from setuptools import setup, find_packages
import os

# https://github.com/readthedocs/readthedocs.org/issues/5512#issuecomment-475073310
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    INSTALL_REQUIRES = []
else:
    INSTALL_REQUIRES = ['obspy']

setup(
      name='waveform_collection',
      packages=find_packages(),
      install_requires=INSTALL_REQUIRES
      package_data={'': ['../avo_json/*.json']},
      )
