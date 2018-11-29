## DupLess (Blast-based Assembly dedupliCation of Heterozygous regions):
## A tool to remove assembly heterozygous duplication based on blast and read coverage.

Most of the currently available assemblers are designed to work on highly inbred, homozygous species and treats differing haplotypes as separate contigs. However, inbreeding is not always an option and attempts to assemble a highly heterozygous species often results in a heavily duplicated assembly. 
For these cases, we created "DupLess”, a tool capable of quickly detecting and removing the duplicated regions issued from heterozygosity.

DupLess workflow is composed of two main steps:
 1. The detection, based on the coverage, of heterozygous regions in the assembly.
 2. The detection of duplicates among the heterozygous regions, based on their sequences similarity (using blast).

---

## Dependencies

- python
- The following python packages: numpy, pandas, subprocess, Bio, matplotlib.pyplot, getopt, multiprocessing.

- **samtools v1.9** or higher
- **bedtools v2.27.1** or higher
- **blastn v2.6.0+** or higher

---

## Input Files:

**Required**

- The assembly in fasta format.
- A bed file with the coverage value for each base of the assembly. You can produce such a file by aligning reads to the assembly and then run "bedtools genomecov" on the resulting bam.

**Optional**

- A bed file contaning the Gaps in the assembly. If provided, they will be represented as grey bars on the graphs.
- If you wish to skip the detection of heterozygous regions based on the coverage, you can directly input a bed file with the regions to consider for duplication. (This file is also produced during DupLess first step)

---

## Usage

	python DupLess.py -t [nb_threads] -w [window_size] -b [coverage.bed] -a [assembly.fasta] -c [expected_coverage] -g [gaps.bed] -i [min_blast_identity] -l [min_blast_length] -o [output_folder]

**Options:**

     -t/--nThreads               The number of threads (default 20)
     -o/--out_folder             The output folder (default is the current directory).

     -b/--bed_cov                The bed file containing the coverage for each position (can be generated with bedtools genomecov).
     -a/--assembly               The assembly corresponding to the bed coverage in fasta format.
     -g/--bed_gaps               A bed file contaning the gaps along the genome. If given, the graphs will contain a grey background where the gaps are.
     
     -w/--window_size            The size of the window. The value of the read coverage will be the median of the values inside each window (default: 1000).
     -c/--expected_cov           The expected read coverage along the genome. The homozygosity / heterozygosity will be determined based on this value.
                                 If no value is given, it will be based on the mode of the coverage distribution (not reliable if high heterozygosity).
     -i/--blast_identity         The minimum percentage of identity between the het region and the blast hit to consider it valid (default: 90, range 0 to 100).
     -l/--blast_length           The minimum length for the blast hit to be considered as valid (default=0).


**Other:**

	-s/--skip_het_detection     Skip the detection of the heterozygous regions. If so, you must provide a bed with the heterozygous regions positions:
                                     python DupLess.py -t [nb_threads] -a [assembly.fasta] -s [het_regions.bed] -i [min_blast_identity] -l [min_blast_length] -o [output_folder]

	-h/--help                   Print the usage and help.

---

## Output

- Two fasta files containing the different versions of the duplications.
- A bed file with the identified heterozygous regions.
- The results of the blast between the heterozygous regions.
- Graphs of the coverage along each sequences of the assembly.
- A histogram of the coverage distribution.
