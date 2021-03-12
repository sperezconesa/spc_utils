"""
Set of python wrappers of the os libray to make things cleaner
and avoid annoing things of os.
"""
import errno
import os
from shutil import copytree


def remove(filename):
    """
    Call remove of os and not give error if file does not exist.
    """
    try:
        os.remove(filename)
    except OSError as e:  # this would be "except OSError, e:" before Python 2.6
        print("Didn't remove anything")
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occurred


def mkdir(filename):
    """
    Call makedirs of os and not give error if path exists
    """
    os.makedirs(filename, exist_ok=True)


def copy_dir(input_path, out_path, verbose=True):
    """
    Copy a directory and if already there do nothing or give a message.
    """
    assert isinstance(verbose, bool), "verbose should be a bool."
    assert os.path.isdir(input_path), "Input  path does not exisist."
    try:
        copytree(input_path, out_path)
    except OSError as e:
        if verbose:
            print("already there")


def natural_sort(l):
    """
    Takes as input a list l of strings and sorts it with natural order.
      Parameters
      ----------
      l: list of strings.
      Returns
      -------
      l sorted
    """
    from re import split

    assert isinstance(l, list), "l is not a list!"
    for i in l:
        assert isinstance(i, str), "List contains non-string elements."
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in split("([0-9]+)", key)]
    return sorted(l, key=alphanum_key)
