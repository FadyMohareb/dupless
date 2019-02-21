#!/usr/bin/python

import sys
import os

def check_return_code(code, cmd):
    if(code != 0):
        print("\nAn error occured during the following command:")
        print(cmd)
        print("Error code:"+str(code))
        sys.exit(code)


def check_file(filename):
    """
    Checks if "filename" exists and is a file.
    Returns:
        True if file exists and is a file.
        False if filename==None or is not a file.
    """
    file_ok = True
    error_mssg = ""
    if(filename == None):
        error_mssg = "Error: file is missing."
        file_ok = False
    else:
        if not os.path.isfile(filename):
            error_mssg = "Error: '"+str(filename)+"' is not a file."
            file_ok = False
    return file_ok, error_mssg
