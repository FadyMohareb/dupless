#!/bin/python

from subprocess import call, Popen, PIPE
from Bio import SeqIO
import os



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



def extract_and_blast(region, het_fasta, output_folder, HetDect_folder):
    """
    Creates commands for blasting a region (format "name:start-stop") against all the other het regions (needs a fasta file with all the het regions "het fasta").
    Returns the commands (strings) to:
        Remove the het regions from the het fasta.
        Extract the het region from the het fasta.
        Create the blast DB from the filtered het fasta.
        Blast the het region against the blast SB.
    """
    region = region.rstrip("\n")
    region_temp_fasta = output_folder+"/temp/"+region+"_temp.fasta"
    region_temp_blastdb = output_folder+"/temp/"+region+"_temp.blastdb"
    # We remove the sequence from the list of het sequences to avoid self hit with blast.
    filtered_fasta = output_folder+"/temp/Het_regions_filtered_"+region+".fasta"

    filter_fasta = "python "+HetDect_folder+"/filter_fasta_by_id.py "+het_fasta+" "+filtered_fasta+" "+region+ " &"
    extract = "samtools faidx "+het_fasta+" "+region+" -o "+region_temp_fasta+" &"
    makeBlastDB = "makeblastdb -in "+filtered_fasta+" -input_type fasta -out "+region_temp_blastdb+" -dbtype nucl &"
    megablast = "blastn -task megablast -db "+region_temp_blastdb+" -query "+region_temp_fasta+" -outfmt 6 -out "+output_folder+"/individual_blasts/"+region+".tab -max_target_seqs 1 -max_hsps 1 &"

    return filter_fasta, extract, makeBlastDB, megablast



def multi_processes(list_cmds):
    """
    Takes a list of commands and launches them in parallel.
    """
    procs = [ Popen(i, shell=True, stdout=PIPE) for i in list_cmds ]
    for p in procs:
        p.communicate()



def detect_dupl_regions(assembly_name, het_bed, output_folder, nbThreads, HetDect_folder):
    """
    Main script to launch the comparison between the het regions.
    Writes the results of Blast in output_folder"/All_Blasts.tab"
    """
    # Extract the heterozygous regions from the assembly and create a fasta file with it.
    het_fasta_name = output_folder+"/assembly_HET_ONLY.fa"
    process = Popen(["bedtools", "getfasta", "-fi", assembly_name, "-bed", het_bed, "-name", "-fo", het_fasta_name], stdout=PIPE)
    process.wait()

    process = Popen(["samtools", "faidx", het_fasta_name], stdout=PIPE)
    process.wait()

    print("Blasting each heterozygous regions against the others, this could take a while....")
    filter_cmds = list()
    extract_cmds = list()       # Extract the region
    makeBlastDB_cmds = list()   # Create a copy of the heterozygous regions fasta without the region and create a database
    megablast_cmds = list()     # Blast the region against the database
    remove_cmds = list()        # Remove the filtered copy
    n_process = 0

    # We do each step by batch to parallelise, all the extractions, then all the makeBlastDB, then all the megablasts
    het_fasta = SeqIO.parse(het_fasta_name, "fasta")
    for seq_record in het_fasta:
        filter_fasta, extract, makeBlastDB, megablast = extract_and_blast(seq_record.name, het_fasta_name, output_folder, HetDect_folder)
        filter_cmds.append(filter_fasta)
        extract_cmds.append(extract)
        makeBlastDB_cmds.append(makeBlastDB)
        megablast_cmds.append(megablast)
        n_process = n_process + 1

        if n_process >= int(nbThreads):
            multi_processes(filter_cmds)
            multi_processes(extract_cmds)
            multi_processes(makeBlastDB_cmds)
            multi_processes(megablast_cmds)

            n_process = 0
            filter_cmds = list()
            extract_cmds = list()
            makeBlastDB_cmds = list()
            megablast_cmds = list()
            process = Popen("rm "+output_folder+"/temp/*", shell=True, stdout=PIPE)
            process.wait()

    multi_processes(filter_cmds)
    multi_processes(extract_cmds)
    multi_processes(makeBlastDB_cmds)
    multi_processes(megablast_cmds)
    multi_processes(remove_cmds)
    print("Blast done !")

    # Concatenate the blast results and filter them by identity% and length.
    concat_blasts_name = output_folder+"/All_Blasts.tab"
    print("Concatenating the blast results to "+concat_blasts_name+" ...")
    with open(concat_blasts_name, "w") as concat_blasts:
        process = Popen(["find", output_folder+"/individual_blasts/", "-maxdepth", "1", "-type", "f", "-exec", "cat", "{}", "+"], stdout=concat_blasts)
        process.wait()
   print("Bed files concatenated to: "+concat_blasts_name)
