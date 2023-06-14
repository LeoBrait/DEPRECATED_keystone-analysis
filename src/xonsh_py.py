#!/usr/bin/env python3

import os
import sys
import shutil
from distutils import dir_util

# defining a python class to redirect script's output to both logfile and stdout

class Logger:
    """
    The code is used for logging the output of the program. This is useful for debugging purposes. The code is used to log the output of the program to a file called log.txt.
    """
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("log.txt", "w")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    def flush(self):
        pass

def existOldData(path=None):
    # This function checks if the old data exists
    # Check if the path is a string
    if not isinstance(path,str):
        raise Exception("Wrong path data type!")
    # Return true if the path exists
    return os.path.exists(path)


def lsgrep(path, wild_regex, with_path=True):
    """
    This code searches for files in the given path that match the given wildcards
    It returns a list of the paths to the matching files
    """
    base = path+'/' if with_path else ''
    ret = [base+f for f in os.listdir(path) if any(map(lambda wild: wild in f, wild_regex))]
    return ret

def mkdir_p(path):
    """
    This code creates a directory if it does not exist
    """
    for p in path:
        os.makedirs(p,mode=0o755,exist_ok=True)

def cp(path1,path2):
    """
    This code copies a file from path1 to path2
    """
    for path in path1:
        shutil.copy2(path,path2)

def cpr(path1,path2):
    """
    This code copies a directory from path1 to path2
    """
    for path in path1:
        dir_util.copy_tree(path, path2)

def rm(path):
    """
    This code removes a file
    """
    for p in path:
        os.remove(p)

def rmr(path):
    """
    This code removes a directory
    """
    for p in path:
        shutil.rmtree(p,ignore_errors=True)

def mv(path1,path2):
    """
    This code moves a file from path1 to path2
    """
    for path in path1:
        shutil.move(path,path2)

def cd(path):
    """
    This code changes the current working directory
    """
    os.chdir(path)

def pwd():
    """
    This code returns the current working directory
    """
    return os.getcwd()

def echo(path,text):
    """
    This code writes text to a file
    """
    with open(path,'w') as f:
        print(text,file=f)

def sexec(cmd):
    """
    This code executes a shell command
    """
    return os.system(cmd)

def cat(path):
    """
    This code reads a file
    """
    with open(path,'r') as f:
        return f.read()
