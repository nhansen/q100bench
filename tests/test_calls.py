import pytest
import subprocess
import sys
from q100bench import bench
from q100bench import seqparse
from q100bench import bedtoolslib

def test_configs():
    args = bench.parse_arguments(['-c', 'tests/testconfig.txt', '-b', 'blah', '-r', 'blah', '-q', 'blah', '-p', 'blah'])
    configvals = bench.read_config_data(args)
    assert configvals['resourcedir'] == 'resources'

def test_genstats():
    args = bench.parse_arguments(['-c', 'tests/testconfig.txt', '-b', 'tests/testbam.bam', '-r', 'tests/testref.fa', '-q', 'tests/testquery.fa', '-p', 'tests/testrun'])
    configvals = bench.read_config_data(args)

