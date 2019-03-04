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


def make_fasta_one_line(fasta_input, fasta_oneLine):
    """
    Transform to single line fasta to avoid empty lines after sed step
    Indeed if a region is longer than the fasta wrapping (usually 80 caracters), then the fasta will contain empty lines.
    awk from: https://stackoverflow.com/questions/15857088/remove-line-breaks-in-a-fasta-file, to avoid error with biopython: "fasta-2line"
    """
    with open(fasta_oneLine, "w") as fasta_oneLine_handle:
        cmd_oneLine = "awk \'/^>/{print s? s\"\\n\"$0:$0;s=\"\";next}{s=s sprintf(\"%s\",$0)}END{if(s)print s}\' "+fasta_input
        print("\t"+cmd_oneLine)
        try:
            pr = subprocess.Popen(cmd_oneLine, shell=True, stdout=fasta_oneLine_handle)
            pr.communicate()
            check_return_code(pr.returncode, cmd_oneLine)
        except:
            print("Error for: " + cmd_oneLine)
            print(sys.exc_info()[0])
            sys.exit()


def remove_file(filename):
    """
    Try to remove a file.
    Do not exit if fails.
    """
    try:
        pr = subprocess.Popen(["rm", filename], shell=False)
        pr.communicate()
        #ud.check_return_code(pr.returncode, "rm "+filename)
    except:
        print("Error for: rm " + filename)
        print(sys.exc_info()[0])


def empty_folder(folder):
    """
    Removes all the files in a folder.
    Used to empty the temp file during the blast step.
    """
    cmd = "rm "+folder+"/*"
    try:
        # The shell=True needed here because of the "*" (regex do not work with shell=False)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.communicate()
    except:
        print("Error for: " + cmd)
        print(sys.exc_info()[0])
        sys.exit()
