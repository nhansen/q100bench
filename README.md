# q100bench

The q100bench python package analyses a user-supplied alignment of a test assembly to a benchmark, and prints general statistics, BED-formatted regions regarding the alignments, and PDF-formatted plots.

The program was written by Nancy Fisher Hansen, a staff scientist in the Genome Informatics Section at the National Human Genome Research Institute (NHGRI). Nancy can be reached at nhansen@mail.nih.gov.

## Install

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

### Dependencies

This program uses R's Rscript command with [Bioconductor](https://www.bioconductor.org/) to create plots, and [bedtools](https://bedtools.readthedocs.io/en/latest/) to compare and merge intervals.

In addition, the program requires various BED-formatted files with data about the benchmark assembly. For the Q100 assembly hg002v1.0.1, a tarball of these files is available on [AWS](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/benchmark/resources/hg002v1.0.1.resources.tar.gz). Once downloaded, the tarball should be unpacked and the locations of files should be included in the config file passed to the program.

All other dependencies are installed by the pip installer with the commands in the previous section. Feel free to post installation issues to the issues section of this github repository.


