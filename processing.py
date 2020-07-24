#!/usr/bin/python

# Part of DupLess dealing with the removal of duplications from the assembly

import subprocess
import sys

import utils_dupless as ud

def remove_duplications_assembly(outname, assembly_name, bedname, output_folder):
    """From an assembly and a bed of regions to remove, create a fasta with duplications removed.
    The regions in the bed are masked as "$" with "bedtools maskfasta", then sed is used to replace "$" by "".
    Also use "awk" to have each sequence on one line (to avoid creating empty lines with sed).
    """
    # Temp files used by sed and awk, removed at the end of the function
    fasta_masked = output_folder+"/deduplicated/temp_masked.fasta"
    fasta_masked_oneLine = output_folder+"/deduplicated/temp_masked_oneLine.fasta"

    # Mask the duplicate regions with "$"
    cmd_mask = ["bedtools", "maskfasta", "-fi", assembly_name, "-fo", fasta_masked, "-bed", bedname, "-mc", "$"]
    print("\t"+ " ".join(cmd_mask))
    try:
        pr = subprocess.Popen(cmd_mask, shell=False, stdout=subprocess.PIPE)
        pr.communicate()
        ud.check_return_code(pr.returncode, " ".join(cmd_mask))
    except Exception as e:
        print("Error for: " + " ".join(cmd_mask))
        print("Exception:"+str(e))
        sys.exit()

    # Make the fasta as one sequence per line, needed for sed below:
    # else we have "empty/truncated" lines where we remove sequences longer than fasta line length.
    ud.make_fasta_one_line(fasta_masked, fasta_masked_oneLine)

    # Replace the "$" which marks the duplicate regions by an empty string (to remove them)
    # The backup extension is required for Mac OS ("-i'backup' -e" should work on both GNU and Mac)
    cmd_sed = "sed -i'.dupless_sed_backup' -e 's/\$//g' "+fasta_masked_oneLine
    print("\t"+cmd_sed)
    try:
        pr = subprocess.Popen(cmd_sed, shell=True, stdout=subprocess.PIPE)
        pr.communicate()
        ud.check_return_code(pr.returncode, cmd_sed)
    except Exception as e:
        print("Error for: " + cmd_sed)
        print("Exception:"+str(e))
        sys.exit()
    
    # Cleaning temp files:
    # move the fasta_masked_oneLine to the deduplicated.fasta
    mv_cmd = ["mv", fasta_masked_oneLine, outname]
    try:
        pr = subprocess.Popen(mv_cmd, shell=False)
        pr.communicate()
        ud.check_return_code(pr.returncode, " ".join(mv_cmd))
    except Exception as e:
        print("Error for: "+ " ".join(mv_cmd))
        print("Exception:"+str(e))
        sys.exit()
    # remove the temp file fasta_masked
    ud.remove_file(fasta_masked)
    # remove the backup made by sed
    ud.remove_file(fasta_masked_oneLine+".dupless_sed_backup")


def generate_discarded_fasta(assembly_name, region_file, discarded_fasta):
    """Write the sequences that have been removed to a discarded.fasta file.
    The bed file should correspond to "discarded.region"
    """
    with open(discarded_fasta, "w") as discarded_handle:
        # Extract the duplicated regions
        cmd_faidx = ["samtools", "faidx", assembly_name, "-r", region_file]
        print("\t"+ " ".join(cmd_faidx))
        try:
            pr = subprocess.Popen(cmd_faidx, shell=False, stdout=discarded_handle)
            pr.communicate()
            ud.check_return_code(pr.returncode, " ".join(cmd_faidx))
        except Exception as e:
            print("Error for: " + " ".join(cmd_faidx))
            print("Exception:"+str(e))
            sys.exit()
