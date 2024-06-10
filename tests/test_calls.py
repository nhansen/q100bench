import pytest
import subprocess
import sys
import os
import pysam
from q100bench import bench
from q100bench import output
from q100bench import seqparse
from q100bench import bedtoolslib

def test_configs():
    args = bench.parse_arguments(['-c', 'tests/testconfig.txt', '-b', 'blah', '-r', 'blah', '-q', 'blah', '-p', 'blah'])
    configvals = bench.read_config_data(args)
    assert configvals['resourcedir'] == 'resources'

def test_checkforprogs():
    bench.check_for_bedtools()
    no_rscript = bench.check_for_R()

def test_createoutputdir():
    args = bench.parse_arguments(['-c', 'tests/testconfig.txt', '-b', 'tests/testsort.bam', '-r', 'tests/testbenchmark.fasta.gz', '-q', 'tests/testassembly.fasta.gz', '-p', 'tests/testrun'])
    configvals = bench.read_config_data(args)
    outputdir = output.create_output_directory(args.prefix)
    assert os.path.isdir(outputdir)

def test_readfasta():
    args = bench.parse_arguments(['-c', 'tests/testconfig.txt', '-b', 'tests/testsort.bam', '-r', 'tests/testbenchmark.fasta.gz', '-q', 'tests/testassembly.fasta.gz', '-p', 'tests/testrun'])
    configvals = bench.read_config_data(args)

    refobj = pysam.FastaFile(args.reffasta)
    queryobj = pysam.FastaFile(args.queryfasta)
    outputfiles = output.name_output_files(args, args.prefix)
    seqparse.write_genome_bedfiles(queryobj, refobj, args, configvals, outputfiles, {})
