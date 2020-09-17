'''
Set of python wrappers of the os libray to make things cleaner
and avoid annoing things of os.
'''
import os, errno

def remove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        print("Didn't remove anything")
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred
def mkdir(filename):
    os.makedirs(filename, exist_ok=True)
