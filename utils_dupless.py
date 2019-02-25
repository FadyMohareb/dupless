#!/usr/bin/python

import sys
import os
import subprocess

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


def index_fasta_file(fasta):
    """
    Index a fasta file with samtools faidx.
    Used for extraction of heterozygous regions.
    """
    cmd = ["samtools", "faidx", fasta]
    print("Generating the index file for the fasta reference:\n\t"+" ".join(cmd))
    try:
        pr = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE)
        pr.communicate()
        check_return_code(pr.returncode, " ".join(cmd))
    except:
        print("Error for: " + " ".join(cmd))
        print(sys.exc_info()[0])
        sys.exit()


def check_old_bedtools_version():
    """
    bedtools before v2.27 adds extra data to fasta sequence names when using getfasta
    This function checks if the current version of bedtools is 2.27 or higher.
    Returns:
        False if bedtools  >= 2.27
        True if bedtools < 2.27
    """
    cmd = ["bedtools", "--version"]
    old = True
    try:
        pr = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE)
        out = str(pr.communicate()[0])
        bedtools_version = out.split("v")[-1]
        main = bedtools_version.split(".")[0]
        second = bedtools_version.split(".")[1]
        if(int(main) >2):
            old = False
        elif(int(main)==2 and int(second)>=27):
            old = False
    except:
        print("Error for: " + " ".join(cmd))
        print(sys.exc_info()[0])
        sys.exit()
    return old