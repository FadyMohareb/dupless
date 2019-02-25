#!/usr/bin/python

import subprocess
from multiprocessing import Pool
from Bio import SeqIO
import os
import sys

import utils_dupless as ud


def get_lengths_fasta(fasta):
    """
    Takes a fasta file name.
    Returns:
        A list containing the length(s) of the sequence(s)
    """
    lengths = []
    fasta_handle = SeqIO.parse(fasta, "fasta")

    for record in fasta_handle:
        lengths.append(len(record.seq))

    return lengths



def extract_and_blast(region, het_fasta, output_folder, DupLess_folder):
    """
    Creates commands for blasting a region (format "name:start-stop") against all the other het regions (needs a fasta file with all the het regions "het fasta").
    Returns the commands (strings) to:
        Remove the het regions from the het fasta.
        Extract the het region from the het fasta.
        Create the blast DB from the filtered het fasta.
        Blast the het region against the blast DB.
    """
    region = region.rstrip("\n")
    region_temp_fasta = output_folder+"/temp/"+region+"_temp.fasta"
    region_temp_blastdb = output_folder+"/temp/"+region+"_temp.blastdb"
    # We remove the sequence from the list of het sequences to avoid self hit with blast.
    filtered_fasta = output_folder+"/temp/Het_regions_filtered_"+region+".fasta"

    filter_fasta = "python "+DupLess_folder+"/filter_fasta_by_id.py "+het_fasta+" "+filtered_fasta+" "+region
    extract = "samtools faidx "+het_fasta+" "+region+" -o "+region_temp_fasta
    makeBlastDB = "makeblastdb -in "+filtered_fasta+" -input_type fasta -out "+region_temp_blastdb+" -dbtype nucl"
    megablast = "blastn -task megablast -db "+region_temp_blastdb+" -query "+region_temp_fasta+" -outfmt 6 -out "+output_folder+"/individual_blasts/"+region+".tab -max_target_seqs 1 -max_hsps 1"

    return filter_fasta, extract, makeBlastDB, megablast



def run_cmd(cmd):
    """
    Used by the multiprocessing.Pool to run commands in parallel
    """
    # Here the try and expect are used to resolve a bug in python, when using multiprocessing and trying to ctrl-c:
    # Subprocesses do not return anything and the pool just creates the next worker, so the ctrl-c does not stop everything
    # We need to catch the keyboard exception here AND in the main thread to exit gracefully
    try:
        subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
    except KeyboardInterrupt:
        return


def get_assembly_coordinates_from_blast_results(region_name, blast_start, blast_stop):
    """
    Translate the blast coordinates (which are relative to the het. regions) to scaffold coordinates.
    The blast hits are relative to the het region coordinates: "region_name start stop" (region_name is in the format: "scaffold_start_stop").
    We do this because we want to remove sequences according to the scaffolds coordinates.
    Used by: "convert_region_coord_to_scaffold_coord"
    Returns:
        The scaffold name in the assembly ("scaffold_start_stop" becomes "scaffold").
        The start position of the blast hit in the assembly "start het region + start blast hit"
        The stop position of the blast hit in the assembly "start het region + blast hit length"
        Basically, we just "shift" the blast coordinates by "region start" to get the scaffold coordinates.
    """
    # We split the region by underscores (name is under the format: "scaffold_start_stop")
    name_split = region_name.split("_")
    # If the original scaffold name contained underscores, we put them back with join 
    # The last 2 elements of the split are the start and end positions of the region
    scaffold_name = "_".join(name_split[:-2])
    region_start = int(name_split[-2])
    #region_end = int(name_split[-1])

    # The part to remove from the scaffold is from: "region start + blast start" to "region start + blast length"
    scaffold_blast_start = region_start + int(blast_start)
    scaffold_blast_end = region_start + int(blast_stop)

    return scaffold_name, scaffold_blast_start, scaffold_blast_end



def convert_region_coord_to_scaffold_coord(region_blasts_name, output_folder):
    """
    Blast hits coordinates are relatives to the heterozygous regions
    But we need the scaffolds coordinates to remove from reference fasta
    This function converts the region coordinates to scaffold coordinates
    """
    scaffold_coord_blasts_name = output_folder+"/All_Blasts_scaffolds_coord.tab"

    print("Converting Blast results in region coordinates to scaffold coordinates...")
    with open(region_blasts_name, "r") as blasts_regions:
        with open(scaffold_coord_blasts_name, "w") as blasts_scaffolds:
            for blast in blasts_regions:
                tabs = blast.split("\t")
                    
                # The QUERY should always be "start stop", but if there is a inverted alignment, the SUBJECT can be "stop   start".
                # So we get the start by taking the min and the stop with the max
                query_start = int(tabs[6])
                query_end = int(tabs[7])

                subject_start = min(int(tabs[8]), int(tabs[9]))
                subject_end = max(int(tabs[8]), int(tabs[9]))

                # We need to get the original positions, since "query_start" and "query_end" correspond to the het regions positions and not the scaffold position.
                q_scaffold_name, q_scaffold_start, q_scaffold_end = get_assembly_coordinates_from_blast_results(tabs[0],query_start, query_end)
                s_scaffold_name, s_scaffold_start, s_scaffold_end = get_assembly_coordinates_from_blast_results(tabs[1], subject_start, subject_end)

                # We want to keep the order of start and end for the subject if there is an reverse alignment.
                if(int(tabs[8]) < int(tabs[9])):
                    blasts_scaffolds.write(q_scaffold_name+"\t"+s_scaffold_name+"\t"+tabs[2]+"\t"+tabs[3]+"\t"+tabs[4]+"\t"+tabs[5]+"\t"+str(q_scaffold_start)+"\t"+str(q_scaffold_end)+"\t"+str(s_scaffold_start)+"\t"+str(s_scaffold_end)+"\t"+tabs[10]+"\t"+tabs[11])
                else:
                    blasts_scaffolds.write(q_scaffold_name+"\t"+s_scaffold_name+"\t"+tabs[2]+"\t"+tabs[3]+"\t"+tabs[4]+"\t"+tabs[5]+"\t"+str(q_scaffold_start)+"\t"+str(q_scaffold_end)+"\t"+str(s_scaffold_end)+"\t"+str(s_scaffold_start)+"\t"+tabs[10]+"\t"+tabs[11])
    print("Blast files with scaffold coordinates written to: "+scaffold_coord_blasts_name+"\n")



def extract_heterozygous_regions(assembly, heterozygous_bed, output_folder):
    """
    Creates a fasta file containing all the heterozygous regions from the "heterozygous_bed" file
    The bed file comes from the previous step of DupLess (or given by the user)
    Returns:
        The name of the fasta file with the heterozygous regions
    """
    het_fasta_name = output_folder+"/assembly_HET_ONLY.fa"

    cmd = ["bedtools", "getfasta", "-fi", assembly, "-bed", heterozygous_bed, "-name", "-fo", het_fasta_name]
    try:
        pr = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE)
        pr.communicate()
        ud.check_return_code(pr.returncode, " ".join(cmd))
    except:
        print("Error for: " + " ".join(cmd))
        print(sys.exc_info()[0])
        sys.exit()
    return het_fasta_name


def concatenate_blast_results(output_folder):
    """
    Concatenate all the individual blasts into one file
    Return:
        The name of the file with concatenated blast results
    """
    region_blasts_name = output_folder+"/All_Blasts_region_coord.tab"
  
    print("Concatenating the blast results to "+region_blasts_name+" ...")
    with open(region_blasts_name, "w") as blasts_regions:
        cmd = ["find", output_folder+"/individual_blasts/", "-maxdepth", "1", "-type", "f", "-exec", "cat", "{}", "+"]
        try:
            pr = subprocess.Popen(cmd, shell=False, stdout=blasts_regions)
            pr.communicate()
            ud.check_return_code(pr.returncode, " ".join(cmd))
        except:
            print("Error for: " + " ".join(cmd))
            print(sys.exc_info()[0])
            sys.exit() 

    return region_blasts_name


def detect_dupl_regions(assembly_name, het_bed, output_folder, nbThreads, DupLess_folder):
    """
    Main script to launch the comparison between the het regions.
    Writes the results of Blast in: "output_folder/individual_blasts/"
    """

    # Extract all the heterozygous regions from the assembly and create a fasta file with it.
    het_fasta_name = extract_heterozygous_regions(assembly_name, het_bed, output_folder)

    # bedtools before v2.27(?) adds the region to the fasta header during getfasta
    # Checks if bedtools < 2.27 and remove 
    if(ud.check_old_bedtools_version()):
        cmd = "sed -i s/::.*//g " + het_fasta_name
        print("bedtools is older than v2.27, removing the part behind the '::'\n\t"+cmd)
        try:
            pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            pr.communicate()
            ud.check_return_code(pr.returncode, " ".join(cmd))
        except:
            print("Error for: " + " ".join(cmd))
            print(sys.exc_info()[0])
            sys.exit()

    ud.index_fasta_file(het_fasta_name)

    print("Blasting each heterozygous regions against the others, this could take a while....")
    filter_cmds = list()
    extract_cmds = list()       # Extract the region
    makeBlastDB_cmds = list()   # Create a copy of the heterozygous regions fasta without the region and create a database
    megablast_cmds = list()     # Blast the region against the database
    remove_cmds = list()        # Remove the filtered copy
    n_process = 0

    pl = Pool(nbThreads)

    # We do each step by batch to parallelise, all the extractions, then all the makeBlastDB, then all the megablasts
    het_fasta = SeqIO.parse(het_fasta_name, "fasta")
    for seq_record in het_fasta:
        # Create the list of commands to use to process 1 region
        filter_fasta, extract, makeBlastDB, megablast = extract_and_blast(seq_record.name, het_fasta_name, output_folder, DupLess_folder)
        filter_cmds.append(filter_fasta)
        extract_cmds.append(extract)
        makeBlastDB_cmds.append(makeBlastDB)
        megablast_cmds.append(megablast)
        n_process = n_process + 1

        # If we get more commands than nbThreads specified, we launch them in parallel
        if n_process >= int(nbThreads):
            # Bug in python: if ctrl-c during multiprocessing, the children becomes unjoinable
            # To resolve this, have to add a keyboardInterrupt exception to both the main thread AND the subprocessess
            try:
                pl.map(run_cmd, filter_cmds)
                pl.map(run_cmd, extract_cmds)
                pl.map(run_cmd, makeBlastDB_cmds)
                pl.map(run_cmd, megablast_cmds)
            except KeyboardInterrupt:
                print("Caught KeyboardInterrupt, terminating workers")
                pl.terminate()
                pl.join()
                sys.exit()
            except:
                print("Error during the processing of the contigs:")
                print(sys.exc_info()[0])
                pl.terminate()
                pl.join()
                sys.exit()

            n_process = 0
            filter_cmds = list()
            extract_cmds = list()
            makeBlastDB_cmds = list()
            megablast_cmds = list()

            cmd = "rm "+output_folder+"/temp/*"
            try:
                # The shell=True needed here because of the "*" (regex do not work with shell=False)
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                process.communicate()
            except:
                print("Error for: " + cmd)
                print(sys.exc_info()[0])
                sys.exit()

    # Do the remaining commands (the ones remaining because "n_process <= int(nbThreads)")
    try:
        pl.map(run_cmd, filter_cmds)
        pl.map(run_cmd, extract_cmds)
        pl.map(run_cmd, makeBlastDB_cmds)
        pl.map(run_cmd, megablast_cmds)
        pl.map(run_cmd, remove_cmds)
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        pl.terminate()
        pl.join()
        sys.exit()
    except:
        print("Error during the processing of the contigs:")
        print(sys.exc_info()[0])
        pl.terminate()
        pl.join()
        sys.exit()

    print("Blast done !\n")
    pl.close()
    pl.join()
    
    # Concatenate the blast results and convert blast coordinates from het regions to coordinates on scaffolds
    region_blasts_name = concatenate_blast_results(output_folder)
    print("Blast files concatenated to: "+region_blasts_name)

    # We need to output the scaffold coordinates for the next step
    # The positions from the blast outputs are relative to the het regions
    convert_region_coord_to_scaffold_coord(region_blasts_name, output_folder)