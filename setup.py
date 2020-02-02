from setuptools import setup, find_packages
import os

# https://github.com/readthedocs/readthedocs.org/issues/5512#issuecomment-475073310
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    INSTALL_REQUIRES = []
    INCLUDE_PACKAGE_DATA = True
else:
    INSTALL_REQUIRES = ['obspy']
    INCLUDE_PACKAGE_DATA = False

setup(
      name='waveform_collection',
      packages=find_packages(),
      install_requires=INSTALL_REQUIRES,
      package_data={'': ['../avo_json/*.json']},
      include_package_data=INCLUDE_PACKAGE_DATA
      )
