import sys
import os
import re
import shutil
import pysam
import argparse
import importlib.resources
from pathlib import Path
from q100bench import errors
from q100bench import output
from q100bench import seqparse
from q100bench import alignparse
from q100bench import phasing
from q100bench import stats
from q100bench import plots

#def check_for_bedtools():
    #if shutil.which("bedtools") is None:
        #sys.print("You don\'t seem to have bedtools in your path. Please install bedtools", file=sys.stderr)
        #exit(1)
    #return 0
#
def check_for_R():
    if shutil.which("Rscript") is None:
        print("You don\'t seem to have Rscript in your path. Plots will not be generated", file=sys.stderr)
        return 1
    return 0

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Print read quality statistics from a bam file of the reads aligned to a benchmark assembly"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('-b', '--bam', required=True, default=None, help='bam file of alignments of the reads to the diploid benchmark')
    parser.add_argument('-r', '--reffasta', type=str, required=True, help='(indexed) fasta file for benchmark reference')
    parser.add_argument('--region', type=str, required=False, default='', help='region of benchmark reference to assess (default is to assess the whole reference genome)')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='prefix for output directory name')
    parser.add_argument('-m', '--minalignlength', type=int, required=False, default=5000, help='minimum length of alignment required to be included in alignment statistics and error counts')
    parser.add_argument('--mincontiglength', type=int, required=False, default=500, help='minimum length for contig to be included in contig statistics')
    parser.add_argument('--nomononucs', action='store_true', required=False, help='skip analysis of mononucleotide length accuracy')
    parser.add_argument('--nobaseerrors', action='store_true', required=False, help='skip analysis of base errors within reads')
    parser.add_argument('-e', '--errorfile', type=str, required=False, default='', help='preexisting file of read errors to report and plot stats for')
    parser.add_argument('-R', '--readsetname', type=str, required=False, default="test", help='name of the assembly being tested--should be query sequence in bam file')
    parser.add_argument('-B', '--benchmark', type=str, required=False, default="truth", help='name of the assembly being used as a benchmark--should be the reference sequence in the bam file')
    parser.add_argument('-c', '--config', type=str, required=False, default="benchconfig.txt", help='path to a config file specifying locations of benchmark data files')

    return parser

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    if not args.bam:
        sys.print("Must specify a bam file with --bam", file=sys.stderr)
        exit(1)

    return args

def read_config_data(args)->dict:
    configfile = args.config
    configpath = Path(configfile)

    configvals = {}
    if not configpath.exists():
        print("Using resource locations from default config file")
        template_res = importlib.resources.files("q100bench").joinpath('benchconfig.txt')
        with importlib.resources.as_file(template_res) as configfile:
            with open(configfile, "r") as cr:
                configline = cr.readline()
                while configline:
                    p = re.compile(r'^([^#\s]+):+\s+(\S+)$')
                    match = p.match(configline)
                    if match:
                        key = match.group(1)
                        value = match.group(2)
                        configvals[key] = value
                    configline = cr.readline()
    else:
        print("Using resource locations from " + configfile)
        with open(configfile, "r") as cr:
            configline = cr.readline()
            while configline:
                p = re.compile(r'^([^#\s]+):+\s+(\S+)$')
                match = p.match(configline)
                if match:
                    key = match.group(1)
                    value = match.group(2)
                    configvals[key] = value
                configline = cr.readline()

    return configvals

def main() -> None:

    args = parse_arguments(sys.argv[1:])
    #check_for_bedtools()
    no_rscript = check_for_R()
    outputdir = output.create_output_directory(args)
    benchparams = read_config_data(args)

    hetsites = phasing.read_hetsites(benchparams["hetsitevariants"])

    alignobj = pysam.AlignmentFile(args.bam, "rb")
    refobj = pysam.FastaFile(args.reffasta)

    outputfiles = output.name_read_stats_files(args, outputdir)

    # evaluate mononucleotide runs:
    if not args.nomononucs:
        print("Assessing accuracy of mononucleotide runs")
        print(benchparams["mononucruns"])
        print(outputfiles["mononucstatsfile"])
        mononucstats = errors.assess_mononuc_read_coverage(alignobj, args.region, benchparams["mononucruns"], outputfiles["mononucstatsfile"], hetsites)
        #stats.write_mononuc_stats(mononucstats, outputfiles, benchmark_stats, args)
    
    # evaluate errors within read alignments:
    if not args.nobaseerrors:
        print("Assessing errors within read alignments")
        print(outputfiles["readerrorfile"])
        errorstats = errors.assess_read_align_errors(alignobj, refobj, outputfiles["readerrorfile"], hetsites, args)
        stats.write_read_error_summary(errorstats, outputfiles)
        #plots.plot_read_error_stats(outputfiles, error_stats, args)

    # plot alignment coverage across assembly and genome:
    #if not no_rscript:
    #    print("Creating plots")
    #    if alignobj is not None:
    #        plots.plot_mononuc_accuracy(args.assembly, outputdir, benchparams["resourcedir"])


if __name__ == "__main__":
    main()
