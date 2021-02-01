'''
Set of python wrappers of the os libray to make things cleaner
and avoid annoing things of os.
'''
import os, errno
from shutil import copytree

def remove(filename):
    '''
    Call remove of os and not give error if file does not exist.
    '''
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        print("Didn't remove anything")
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred
def mkdir(filename):
    '''
    Call makedirs of os and not give error if path exists
    '''
    os.makedirs(filename, exist_ok=True)

def copy_dir(input_path, out_path,verbose=True):
    '''
    Copy a directory and if already there do nothing or give a message.
    '''
    assert isinstance(verbose, bool), 'verbose should be a bool.'
    assert os.path.isdir(input_path) ,'Input  path does not exisist.'
    try:
        copytree(input_path, out_path)
    except OSError as e:
        if verbose:
            print('already there')

