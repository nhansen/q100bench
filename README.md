# q100bench

The q100bench python package analyses a user-supplied alignment of a test assembly to a benchmark assembly (preferably from the same sample), and prints general statistics, BED-formatted regions regarding the alignments and discrepancies within them, and PDF-formatted plots.

The program was written by Nancy Fisher Hansen, a staff scientist in the Genome Informatics Section at the National Human Genome Research Institute (NHGRI). Nancy can be reached at nhansen@mail.nih.gov.

## Table of contents
- [Install](#install)
- [Getting Started](#getting-started)
- [Outputs](#outputs)

# Install

## Dependencies

This program uses R's Rscript command with [Bioconductor](https://www.bioconductor.org/) to create plots, and [bedtools](https://bedtools.readthedocs.io/en/latest/) to compare and merge intervals. If the "Rscript" command is not in a user's path, the program will complain and then skip all plotting. If the "bedtools" command isn't in the user's path, the program exits with an error.

In addition, the program requires various BED-formatted files with data about the benchmark assembly. For the Q100 benchmark assembly hg002v1.0.1, a tarball of these files is available on [AWS](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/benchmark/resources/hg002v1.0.1.resources.tar.gz). Once downloaded, the tarball should be unpacked and the locations of files should be included in the config file passed to the program.

All other dependencies are installed by the pip installer with the commands in the next section called "Local Installation". Feel free to post installation issues to the issues section of this github repository.

## Local Installation

Until q100bench is available on PyPi and bioconda, the easiest way to use it is to install it locally. First clone this github repository:
```
git clone https://github.com/nhansen/q100bench
cd q100bench
```

Create a virtual environment for the project:
```
python3 -m venv venv
source venv/bin/activate
```

Finally use python's pip installer to install and test a development copy for yourself to run:
```
python3 -m pip install -e .
pytest
```
# Getting Started

## Evaluating diploid assemblies

At this point, the q100bench program is written to evaluate a BAM-formatted file of alignments of a single haplotype of a diploid assembly to a diploid benchmark genome (e.g., hg002v1.0.1). We recommend aligning each haplotype of the diploid assembly separately to the diploid benchmark using minimap2 with the "-x asm5" preset. To evaluate each haplotype the program usage is

	q100bench -b <assemblyhaplotype_vs_benchmark.bam> -r <benchmark.fasta> -q <assemblyhaplotype.fasta> -p <prefix_for_output> -A <assemblyhaplotype_name> -B <benchmark_name>

For typical assemblies, the program will use about 16Gb of memory and around 15 minutes of CPU time. The command "q100bench --help" will display information on other options available (e.g., to restrict regions of the genome examined, set minimum contig or alignment lengths for processing, etc.).

## Evaluating read sets

To report and plot statistics about discrepancies between a set of sequencing reads and a benchmark diploid genome, the program has a "readbench" command. First, the reads should be aligned to the diploid benchmark assembly with whatever aligner and parameters you feel are most accurate. The usage of the readbench command is

	readbench -b <reads_vs_benchmark.bam> -r <benchmark.fasta> -p <prefix_for_output> -B <benchmark_name> -R <readset_name>

Because it is evaluating more alignments than for an assembly evaluation, the readbench command takes longer to run than the q100bench command. For this reason, it has a "--downsample" option which allows the user to pass a fraction between 0 and 1.0 that will cause read alignments to be randomly downsampled to include only that fraction of the alignments in its accuracy calculations. As with the q100bench command, information about options can be obtained with "readbench --help".


