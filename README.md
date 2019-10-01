waveform_collection
===================

Package containing convenience functions to collect seismic/infrasound waveforms
and metadata from IRIS/WATC/AVO servers or local miniSEED files.

Installation
------------
To use _waveform_collection_, first clone or download this repository. Then
create a new [conda](https://docs.conda.io/projects/conda/en/latest/index.html)
environment with the necessary dependencies (or use a pre-existing one):
```
$ conda env create -n waveform_collection -c conda-forge obspy
```
Then execute
```
$ conda activate waveform_collection
$ cd /path/to/waveform_collection
$ pip install -e .
```
to install the package.

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
