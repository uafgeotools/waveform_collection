#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .common import check_file_exists
import os
from shutil import copyfile
from obspy import read
import subprocess
import tempfile

def read_wfdisc(wfdisc_file):
    '''
    st = read_wfdisc(wfdisc_file)
    
    read_wfdisc takes the path to a CSS wfdisc file and uses the obspy.read(wfdisc_file,'css')
    to read the data. In likely event that the path to the binary data file (mwf) in
    the wfdisc file does not exist, a temporary wfdisc file is used and the path is
    changed (using sed) to match the path to the input wfdisc file.
    
    '''
    
    check_file_exists(wfdisc_file)
    try:
        st = read(wfdisc_file,'css')
    except FileNotFoundError as Err:
        print('fixing mwf file path in wfdisc...')
        mwf_file = '.'.join([os.path.splitext(wfdisc_file)[0],'mwf'])
        #fixed_wfdisc_file = os.path.join(os.sep,'tmp',os.path.basename(wfdisc_file))
        ntf = tempfile.NamedTemporaryFile()
        fixed_wfdisc_file = ntf.name
        copyfile(wfdisc_file,fixed_wfdisc_file)
        
        _fix_mwr_path(fixed_wfdisc_file, os.path.dirname(Err.filename), os.path.dirname(mwf_file))
        
        st = read(fixed_wfdisc_file,'css')
        ntf.close()
    return st

def _fix_mwr_path(wf, match, replacement):
    '''
    fixes the binary data file path in a given wfdisc file using sed
    
    if the match is longer than what needs to be there, pads the sed replacement with spaces;
    if the match is shorter than what needs to be there, pads the sed match with spaces;
    if the match is the same length as what needs to be there, one-for-one sed will do
    '''
    
    diff = len(match) - len(replacement)
    
    if diff > 0:
        # pad replacement with spaces
        replacement = replacement + ' '*diff
    elif diff < 0:
        # pad match with spaces
        match = match + ' '*diff
    elif diff == 0:
        # we're lucky, we don't have to do anything here
        pass
    
    # the maximum length for the Dir field in the wfdisc table is 64 characters
    assert (len(match) < 64) and (len(replacement) < 64)
    
    command = "sed -i -e 's%{}%{}%' {}".format(match,replacement,wf)
    subprocess.call([command], shell=True)