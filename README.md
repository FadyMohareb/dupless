# DupLess v1.0.0 (Duplication Less):

## A tool to remove assembly heterozygous duplication based on blast and read coverage.

Most of the currently available assemblers are designed to work on highly inbred, homozygous species and are treating differing haplotypes as separate contigs. However, inbreeding is not always an option and attempts to assemble a highly heterozygous species often results in a heavily duplicated assembly.
For these cases, we created "DupLess", a tool capable of detecting and removing the duplicated regions issued from heterozygosity.

DupLess workflow is composed of two main steps:

 1. The detection, based on the coverage, of heterozygous regions in the assembly.

 2. The detection of duplicates among the heterozygous regions, based on their sequences similarity (using blast).

---

## Dependencies

- **Python v2 or v3**
- The following python packages: numpy, pandas, biopython, matplotlib.pyplot, getopt, subprocess, multiprocessing, sys, os.

- **samtools v1.9** or higher (/!\ DupLess will not work with version prior to 1.9, as it needs the "-o" parameter)
- **bedtools v2.27.1** or higher (lower versions should also work now, but only v2.26 has been tested)
- **blastn v2.6.0+** or higher
- **awk and sed**


## Installation

- DupLess itself is a collection of python scripts, no installation is needed.

- To install samtools 1.9 (version 1.9 is not yet available from "apt-get install"):
```
	wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
	tar -vxjf samtools-1.9.tar.bz2
	cd samtools-1.9
	make
	export PATH=/path/to/samtools/:$PATH
```
For samtools you may also have to install HTSlib (cf: https://www.biostars.org/p/328831/)

---

## Input Files:

**Required**

- The assembly in fasta format.

- A bed file with the coverage value for each base of the assembly.

**Optional**

- A bed file containing the gaps coordinates in the assembly. If provided, they will be represented as grey bars on the coverage graphs.

- If you wish to skip the detection of heterozygous regions based on the coverage, you can directly input a bed file with the regions to consider for duplication. (This file is also produced during DupLess first step)

**How to generate the coverage bed file**

You need to generate a file with the coverage value at each position (format: "sequence_name   position  coverage"). You can use any pipeline you want to generate this file.

We are using the following pipeline to generate the coverage file (/!\ do not forget the "-d" option for "genomecov" as we need the coverage on every bases):
```
bwa index genome.fa
bwa mem -t 20 genome.fa read_1.fq read_2.fq | samtools view -b -@ 20 -o genome_reads.bam
samtools sort -@ 20 genome_reads.bam > genome_reads.sorted.bam
bedtools genomecov -ibam genome_reads.sorted.bam -d > genome_reads.coverage
```

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

	-n/--no_plot		        Skip the creation of all plots.

	-s/--skip_het_detection     Skip the detection of the heterozygous regions. If so, you must provide a bed with the heterozygous regions positions:
                                     python DupLess.py -s [het_regions.bed] -t [nb_threads] -a [assembly.fasta] -i [min_blast_identity] -l [min_blast_length] -o [output_folder]

	-h/--help                   Print the usage and help and exit.
     -v/--version                Print the version and exit.


## Output

- **Two fasta files containing the different versions of the deduplicated assembly.**
- A bed file with the identified heterozygous regions (useful for exploration of the regions).
- A histogram of the coverage distribution, to help the user decide the expected coverage value (see below).
- Graphs of the coverage along each sequences of the assembly (see below).
- The results of the blast between the heterozygous regions.

![alt text](https://bitbucket.org/MCorentin/hetdect/src/master/exemple_output/Histogram_coverage.png "Histogram of coverage")

![alt text](../exemple_output/Super_scaffold_1.png "Graph of coverage along a sequence")

---

## Tests

To add...

---

## Future work

- Support long reads to realign to duplicated regions and detect misassemblies
- Possibly: add Mummer to check how the regions around the duplications align to each other.
- Add an option to remove whole contigs if the blast hit span > threshold% of the total length.
- Improve speed.
- Flag regions with half the coverage but no blast hits.

## Performances

For a genome of 1Gbp, DupLess requires ~50 Gb of RAM.

For running time, the bottleneck is the pairwise blasting of heterozygous regions, an assembly with a lot of regions will increase the running time significantly.

---