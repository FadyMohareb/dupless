#!/bin/python

from subprocess import call, Popen, PIPE
from Bio import SeqIO
import sys

import utils_dupless as ud

def usage():
    print("python filter_fasta_by_id.py original_fasta new_fasta ID_to_remove")

# ============================ check input ================================
# Only 3 arguments needed, but length of 4 because sys.argv[0] = name of the script
if(len(sys.argv) < 4):
    usage()
    print("3 arguments expected, only "+str(len(sys.argv)-1)+ " found")
    sys.exit(2)

old_fasta_name = sys.argv[1]
new_fasta_name = sys.argv[2]
ID = sys.argv[3]


file_ok, error_mssg = ud.check_file(old_fasta_name)
if not(file_ok):
    print(error_mssg)
    sys.exit(2)

# ======================= Main ==================================
old_fasta = SeqIO.parse(old_fasta_name, "fasta")
ID = ID.rstrip("\n")
ID_found = False

with open(new_fasta_name, "w") as new_fasta:
    for seq_record in old_fasta:
        # If match, do nothing, else write the sequence
        if(str(seq_record.name) != ID):
            SeqIO.write(seq_record, new_fasta, "fasta")
        else:
            ID_found = True

if not ID_found:
    print("Warning: ID '"+ID+"' not found.")
