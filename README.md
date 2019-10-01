waveform_collection
===================

Script containing convenience functions to collect seismic/infrasound waveforms
and metadata from IRIS/WATC/AVO servers or local miniSEED files.

Installation
------------
To use _waveform_collection_, first clone or download this repository. Then the
script can be installed using Python [`setuptools`](https://pypi.org/project/setuptools/),
which allows for installation into your chosen environment by your preferred
method. For example, installing with [pip](https://pypi.org/project/pip/) into a
pre-existing conda environment:
```
$ conda activate my_env
$ cd /path/to/waveform_collection
$ pip install -e .
```
With this method, dependencies (detailed below) are automatically installed if
required.

Dependencies
------------

Python packages:

* [ObsPy](http://docs.obspy.org/)

...and its dependencies, which you don't really have to be concerned about if
you're using [conda](https://docs.conda.io/projects/conda/en/latest/index.html)!

Usage
-----

Access the script's functions with (for example)
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
