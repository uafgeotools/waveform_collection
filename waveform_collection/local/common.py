# -*- coding: utf-8 -*-

import os
import json

def check_file_exists(file):
    '''
    checks a that the specified file is present and readable
    
    Inputs:
        file    [string]    path to file to test
    Outputs:
        exists  [bool]      True if file exists and is readable
    '''
    
    try:
        with open(file,'r'):
            return True
    except Exception as E:
        raise E


def check_file_extension(file, extension):
    '''checks that the specified file has the correct extension'''
    
    if os.path.splitext(file)[-1].lower() == extension.lower():
        return True
    else:
        return False


def load_json_file(file):
    '''loads the specified json file (must have .json extension) using json.load()'''
    
    file = os.path.abspath(file)
    assert check_file_exists(file), 'Error: file {} does not exists or is not readable!'.format(file)
    assert check_file_extension(file, '.json'), 'Error: file {} does not have extension \'{}\'!'.format(file, '.json')
    
    with open(file,'r') as f:
        contents = json.load(f)
        
    return contents