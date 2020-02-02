# -*- coding: utf-8 -*-

import os
import json

def check_file_exists(file):
    """
    Checks that the specified file is present and readable.

    Args:
        file (str): Path to file to test

    Returns:
        bool: `True` if file exists and is readable, `False` otherwise
    """

    try:
        with open(file,'r'):
            return True
    except Exception as E:
        raise E


def check_file_extension(file, extension):
    """
    Checks that the specified file has the correct extension.

    Args:
        file (str): Path to file to test
        extension (str): Extension to test for, e.g. `'.txt'`

    Returns:
        bool: `True` if `file` has extension `extension`, `False` otherwise
    """

    if os.path.splitext(file)[-1].lower() == extension.lower():
        return True
    else:
        return False


def load_json_file(file):
    """
    Loads the specified JSON file (must have ``.json`` extension) using
    :func:`json.load`.

    Args:
        file (str): Path to JSON file

    Returns:
        The contents of the JSON file
    """

    file = os.path.abspath(file)
    assert check_file_exists(file), 'Error: file {} does not exists or is not readable!'.format(file)
    assert check_file_extension(file, '.json'), 'Error: file {} does not have extension \'{}\'!'.format(file, '.json')

    with open(file,'r') as f:
        contents = json.load(f)

    return contents
