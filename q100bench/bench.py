import sys
import os
import re
import shutil
import pysam
import argparse
from pybedtools import BedTool
import importlib.resources
from pathlib import Path
from collections import namedtuple
from q100bench import bedtoolslib
from q100bench import errors
from q100bench import output
from q100bench import seqparse
from q100bench import alignparse
from q100bench import phasing
from q100bench import stats
from q100bench import mummermethods
from q100bench import plots

# create namedtuple for bed intervals:
bedinterval = namedtuple('bedinterval', ['chrom', 'start', 'end', 'name', 'rest']) 


def check_for_bedtools():
    if shutil.which("bedtools") is None:
        sys.print("You don\'t seem to have bedtools in your path. Please install bedtools", file=stderr)
        exit(1)
    return 0

def check_for_R():
    if shutil.which("Rscript") is None:
        print("You don\'t seem to have Rscript in your path. Plots will not be generated", file=stderr)
        return 1
    return 0

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Print assembly statistics from bam files of the assembly aligned to a benchmark assembly"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('-b', '--bam', required=False, default=None, help='bam file of alignments of the test (haploid) assembly to the diploid benchmark')
    parser.add_argument('--paf', required=False, default=None, help='paf-formatted file of alignments of the test (haploid) assembly to the diploid benchmark')
    parser.add_argument('-r', '--reffasta', type=str, required=True, help='(indexed) fasta file for benchmark reference')
    parser.add_argument('-q', '--queryfasta', type=str, required=True, help='(indexed) fasta file for test assembly')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='prefix for output directory name')
    parser.add_argument('-m', '--minalignlength', type=int, required=False, default=5000, help='minimum length of alignment required to be included in alignment statistics and error counts')
    parser.add_argument('--mincontiglength', type=int, required=False, default=500, help='minimum length for contig to be included in contig statistics')
    parser.add_argument('--minns', type=int, required=False, default=10, help='minimum number of consecutive Ns required to break scaffolds into contigs')
    parser.add_argument('-n', '--n_bedfile', type=str, required=False, default=None, help='pre-existing bedfile of locations of N-stretches splitting scaffolds into contigs')
    parser.add_argument('--variantfile', type=str, required=False, default=None, help='pre-existing file of variant locations in assembly compared to benchmark')
    parser.add_argument('-A', '--assembly', type=str, required=False, default="test", help='name of the assembly being tested--should be query sequence in bam file')
    parser.add_argument('-B', '--benchmark', type=str, required=False, default="truth", help='name of the assembly being used as a benchmark--should be the reference sequence in the bam file')
    parser.add_argument('-c', '--config', type=str, required=False, default="benchconfig.txt", help='path to a config file specifying locations of benchmark data files')

    return parser

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    if not args.bam and not args.paf:
        sys.print("Must specify either a bam file with --bam or a paf file with --paf", file=sys.stderr)
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
    check_for_bedtools()
    no_rscript = check_for_R()
    outputdir = output.create_output_directory(args)
    benchparams = read_config_data(args)

    hetsites = phasing.read_hetsites(benchparams["hetsitevariants"])
    hetarrays = phasing.sort_chrom_hetsite_arrays(hetsites)

    alignobj = None
    pafaligns = None
    if args.bam:
        alignobj = pysam.AlignmentFile(args.bam, "rb")
        #bamaligns = alignparse.read_bam_aligns(alignobj, args.minalignlength)
        #qlis_aligns = mummermethods.filter_aligns(bamaligns, "query")
        #rlis_aligns = mummermethods.filter_aligns(bamaligns, "target")
    else:
        pafaligns = alignparse.read_paf_aligns(args.paf, args.minalignlength)
        #qlis_aligns = mummermethods.filter_aligns(pafaligns, "query")
        #rlis_aligns = mummermethods.filter_aligns(pafaligns, "target")

    refobj = pysam.FastaFile(args.reffasta)
    queryobj = pysam.FastaFile(args.queryfasta)

    outputfiles = output.name_output_files(args, outputdir)

    print("Writing bed files for scaffold spans, lengths, N stretches, and contigs (ignoring stretches of less than " + str(args.minns) + " Ns)")
    testgenomebed = seqparse.write_genome_bedfile(queryobj, args, genome_bedfile=outputfiles["testgenomebed"])
    [testcontigbed, testgapbed] = seqparse.find_all_ns(queryobj, args, n_bedfile=outputfiles["testnbed"], atgc_bedfile=outputfiles["testnonnbed"])
    
    print("Writing general statistics about " + args.assembly + " Assembly")
    benchmark_stats = stats.write_general_assembly_stats(outputfiles["generalstatsfile"], refobj, queryobj, testcontigbed, testgapbed, args)

    print("Writing bed files of regions covered by alignments of " + args.assembly + " to " + args.benchmark)
    [refcoveredbed, querycoveredbed, variants, hetsitealleles] = alignparse.write_bedfiles(alignobj, pafaligns, refobj, queryobj, hetarrays, outputfiles["testmatcovered"], outputfiles["testpatcovered"], outputfiles["truthcovered"], outputfiles["variantbed"], outputfiles["coveredhetsitealleles"], args)

    # create merged unique outputfiles:
    [mergedtruthcoveredbed, outputfiles["mergedtruthcovered"]] = bedtoolslib.mergebed(outputfiles["truthcovered"])
    [mergedtestmatcoveredbed, outputfiles["mergedtestmatcovered"]] = bedtoolslib.mergebed(outputfiles["testmatcovered"])
    [mergedtestpatcoveredbed, outputfiles["mergedtestpatcovered"]] = bedtoolslib.mergebed(outputfiles["testpatcovered"])

    print("Writing primary alignment statistics about " + args.assembly + " assembly")
    benchmark_stats = stats.write_merged_aligned_stats(refobj, queryobj, mergedtruthcoveredbed, mergedtestmatcoveredbed, mergedtestpatcoveredbed, outputfiles, benchmark_stats, args)
    #benchmark_stats = stats.write_aligned_stats(refobj, queryobj, refcoveredbed, querycoveredbed, outputfiles, benchmark_stats, args)

    if alignobj is not None:
        # classify variant errors as phasing or novel errors:
        print("Writing phase switch statistics")
        stats.write_het_stats(outputfiles, benchmark_stats, args)
        print("Determining whether errors are switched haplotype or novel")
        benchmark_stats = errors.classify_errors(refobj, queryobj, variants, hetsites, outputfiles, benchparams, benchmark_stats, args)
        stats.write_qv_stats(benchmark_stats, outputfiles, args)

        # measure het phasing across chromosomes:
        #print("Evaluating phasing of heterozygous sites across chromosomes")
        #[coveredhetbed, coveredhetfile] = bedtoolslib.intersectbed(benchparams["hetsitevariants"], outputfiles["truthcovered"], outputfile=outputfiles["coveredhetsitealleles"], writefirst=True)
        #phasing.assess_het_sites(alignobj, refobj, queryobj, coveredhetbed, args)

        # evaluate mononucleotide runs:
        print("Assessing accuracy of mononucleotide runs")
        bedtoolslib.intersectbed(benchparams["mononucruns"], outputfiles["mergedtruthcovered"], outputfile=outputfiles["coveredmononucsfile"], writefirst=True)
        mononucswithvariantsbedfile = bedtoolslib.intersectbed(outputfiles["coveredmononucsfile"], outputfiles["bencherrortypebed"], outputfiles["mononucswithvariantsfile"], outerjoin=True, writeboth=True)
        mononucstats = errors.gather_mononuc_stats(outputfiles["mononucswithvariantsfile"], outputfiles["mononucstatsfile"])
        stats.write_mononuc_stats(mononucstats, outputfiles, benchmark_stats, args)
    
    # plot alignment coverage across assembly and genome:
    if not no_rscript:
        print("Creating plots")
        plots.plot_benchmark_align_coverage(args.assembly, args.benchmark, outputdir, benchparams["resourcedir"])
        plots.plot_testassembly_align_coverage(args.assembly, outputdir, benchparams["resourcedir"])
        if alignobj is not None:
            plots.plot_mononuc_accuracy(args.assembly, outputdir, benchparams["resourcedir"])


if __name__ == "__main__":
    main()
