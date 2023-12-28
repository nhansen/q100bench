import sys
import os
import re
import shutil
import pysam
import argparse
from pybedtools import BedTool
from pathlib import Path
from collections import namedtuple

#sys.path.insert(0, '/data/Phillippy/projects/HG002_diploid/benchmarking/software')

from q100bench import bedtoolslib
from q100bench import errors
from q100bench import output
from q100bench import seqparse
from q100bench import alignparse
from q100bench import stats
from q100bench import plots

# create namedtuple for bed intervals:
bedinterval = namedtuple('bedinterval', ['chrom', 'start', 'end', 'name', 'rest']) 


def check_for_bedtools():
    if shutil.which("bedtools") is None:
        print("You don\'t seem to have bedtools in your path. Please install bedtools")
        exit(1)
    return 0

def check_for_R():
    if shutil.which("Rscript") is None:
        print("You don\'t seem to have Rscript in your path. Plots will not be generated")
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
    parser.add_argument('-b', '--bam', required=True, help='bam file of alignments to search for perfectly aligned segments')
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

def parse_arguments():
    parser = init_argparse()
    args = parser.parse_args()

    return args

def read_config_data(args)->dict:
    configfile = args.config
    configpath = Path(configfile)

    configvals = {}
    with configpath.open("r") as cr:
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

    args = parse_arguments()
    check_for_bedtools()
    no_rscript = check_for_R()
    outputdir = output.create_output_directory(args)
    benchparams = read_config_data(args)

    bamfile = args.bam
    reffasta = args.reffasta
    queryfasta = args.queryfasta
    alignobj = pysam.AlignmentFile(bamfile, "rb")
    refobj = pysam.FastaFile(reffasta)
    queryobj = pysam.FastaFile(queryfasta)

    outputfiledict = output.name_output_files(args, outputdir)

    print("Writing bed files for scaffold lengths, N stretches, and contigs (ignoring stretches of less than " + str(args.minns) + " Ns)")
    genomeintervals = seqparse.write_genome_bedfile(queryobj, args, genome_bedfile=outputfiledict["testgenomebed"])
    [contiglist, gaplist] = seqparse.find_all_ns(queryobj, args, n_bedfile=outputfiledict["testnbed"], atgc_bedfile=outputfiledict["testnonnbed"])
    
    print("Writing general statistics about " + args.assembly + " Assembly")
    benchmark_stats = stats.write_general_assembly_stats(outputfiledict["generalstatsfile"], refobj, queryobj, contiglist, gaplist, args)

    print("Writing " + outputfiledict["truthcovered"] + ", " + outputfiledict["testpatcovered"] + " and " + outputfiledict["testmatcovered"])
    [refcovered, querycovered, variants] = alignparse.write_bedfiles(alignobj, refobj, queryobj, outputfiledict["testmatcovered"], outputfiledict["testpatcovered"], outputfiledict["truthcovered"], outputfiledict["variantbed"], args)

    # create merged unique outputfiles:
    [mergedtruthintervals, outputfiledict["mergedtruthcovered"]] = bedtoolslib.mergebed(outputfiledict["truthcovered"])
    [mergedtestmatintervals, outputfiledict["mergedtestmatcovered"]] = bedtoolslib.mergebed(outputfiledict["testmatcovered"])
    [mergedtestpatintervals, outputfiledict["mergedtestpatcovered"]] = bedtoolslib.mergebed(outputfiledict["testpatcovered"])

    print("Writing primary alignment statistics about " + args.assembly + " assembly")
    benchmark_stats = stats.write_aligned_stats(refobj, queryobj, refcovered, mergedtruthintervals, mergedtestmatintervals, mergedtestpatintervals, outputfiledict, benchmark_stats, args)

    # classify variant errors as phasing or novel errors:
    print("Determining whether errors are switched haplotype or novel")
    variantfile = errors.classify_errors(refobj, queryobj, refcovered, querycovered, variants, outputfiledict, benchparams, args)
    stats.write_qv_stats(benchmark_stats, variantfile, outputfiledict, args)

    # measure het phasing across chromosomes:
    print("Evaluating phasing of heterozygous sites across chromosomes")

    # evaluate mononucleotide runs:
    print("Assessing accuracy of mononucleotide runs")
    bedtoolslib.intersectbed(benchparams["mononucruns"], outputfiledict["mergedtruthcovered"], outputfile=outputfiledict["coveredmononucsfile"], writefirst=True)
    mononucswithvariantsbedfile = bedtoolslib.intersectbed(outputfiledict["coveredmononucsfile"], outputfiledict["bencherrortypebed"], outputfiledict["mononucswithvariantsfile"], outerjoin=True, writeboth=True)
    mononucstats = errors.gather_mononuc_stats(outputfiledict["mononucswithvariantsfile"], outputfiledict["mononucstatsfile"])
    stats.write_mononuc_stats(mononucstats, outputfiledict, benchmark_stats, args)
    
    # plot alignment coverage across assembly and genome:
    if not no_rscript:
        print("Creating plots")
        plots.plot_benchmark_align_coverage(args.assembly, outputdir, benchparams["resourcedir"])
        plots.plot_testassembly_align_coverage(args.assembly, outputdir, benchparams["resourcedir"])
        plots.plot_mononuc_accuracy(args.assembly, outputdir, benchparams["resourcedir"])


if __name__ == "__main__":
    main()
