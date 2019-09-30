waveform_collection
===================

Script containing convenience functions to collect seismic/infrasound waveforms
and metadata from IRIS/WATC/AVO servers or local miniSEED files.


Installation
------------
This package can be installed using python [setuptools](https://pypi.org/project/setuptools/),
which allows installation into your environment by your preferred method. For example,
[pip](https://pypi.org/project/pip/) can be used.

```
cd waveform_collection
pip install .
```


Dependencies
------------

Python packages:

* [ObsPy](http://docs.obspy.org/)

...and its dependencies, which you don't really have to be concerned about if
you're using [conda](https://docs.conda.io/projects/conda/en/latest/index.html)!

It's recommended that you create a new conda environment to use with this
repository:
```
conda create -n waveform_collection -c conda-forge obspy
```

Usage
-----

To use _waveform_collection_, clone or download this repository and add it to
your `PYTHONPATH`, e.g. in a script where you'd like to use
_waveform_collection_:
```python
import sys
sys.path.append('/path/to/waveform_collection')
```
Then you can access package functions with (for example)
```python
from waveform_collection import gather_waveforms
```
and so on. Currently, documentation only exists in function docstrings.

Authors
-------

(_Alphabetical order by last name._)

David Fee  
Liam Toney
Andrew Winkelman
