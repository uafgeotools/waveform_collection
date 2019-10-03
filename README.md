waveform_collection
===================

Package containing convenience functions to collect seismic/infrasound waveforms
and metadata from IRIS/WATC/AVO servers or local files (miniSEED, etc.).

Installation
------------

It's recommended that you install this package into a new or pre-existing
[conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment.
(If you choose the latter option, ensure that your environment contains all of
the packages listed in the [Dependencies](#dependencies) section.)

To create a new conda environment for use with this and other _uafgeotools_
packages, execute the following terminal command:
```
$ conda env create -n uafinfra -c conda-forge obspy
```
This creates a new environment called `uafinfra` with ObsPy and its dependencies
installed.

To install _waveform_collection_, execute the following terminal commands:
```
$ conda activate uafinfra  # Or your pre-existing env
$ git clone https://github.com/uafgeotools/waveform_collection.git
$ cd waveform_collection
$ pip install -e .
```
The final command installs the package in "editable" mode, which means that you
can update it with a simple `git pull` in your local repository. This install
command only needs to be run once.

Dependencies
------------

Python packages:

* [ObsPy](http://docs.obspy.org/)

...and its dependencies, which you don't really have to be concerned about if
you're using conda!

Usage
-----

Access the package's functions with (for example)
```python
from waveform_collection import gather_waveforms
```
and so on. Currently, documentation only exists in function docstrings. For a
usage example, see [`example.py`](example.py).

Authors
-------

(_Alphabetical order by last name._)

David Fee  
Liam Toney  
Andrew Winkelman
