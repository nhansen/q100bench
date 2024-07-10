import sys
import os
import re
import shutil
import pysam
import argparse
import logging
from pybedtools import BedTool
import importlib.resources
from pathlib import Path
from collections import namedtuple
from q100bench import bedtoolslib
from q100bench import errors
from q100bench import output
from q100bench import seqparse
from q100bench import alignparse
from q100bench import structvar
from q100bench import phasing
from q100bench import stats
from q100bench import mummermethods
from q100bench import plots

# create namedtuple for bed intervals:
varianttuple = namedtuple('varianttuple', ['chrom', 'start', 'end', 'name', 'vartype', 'excluded']) 

logger = logging.getLogger(__name__)

def check_for_bedtools():
    if shutil.which("bedtools") is None:
        logger.critical("You don\'t seem to have bedtools in your path. Please install bedtools")
        exit(1)
    return 0

def check_for_R():
    if shutil.which("Rscript") is None:
        logger.warning("You don\'t seem to have Rscript in your path. Plots will not be generated")
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
    parser.add_argument('--maxclusterdistance', type=int, required=False, default=10000, help='maximum distance within a cluster of alignments')
    parser.add_argument('--includefile', type=str, required=False, default=None, help='bed file of benchmark locations to include in the evaluation (stretches of 10 or more Ns and regions excluded in the exclude file will still not be considered)')
    parser.add_argument('--excludefile', type=str, required=False, default=None, help='bed file of benchmark locations to exclude from consideration (in addition to stretches of 10 or more Ns and regions in the exclude file specified in the config file)')
    parser.add_argument('--vcf', action='store_true', required=False, default=False, help='write differences from benchmark in VCF format')
    parser.add_argument('-n', '--n_bedfile', type=str, required=False, default=None, help='pre-existing bedfile of locations of N-stretches splitting scaffolds into contigs')
    parser.add_argument('--variantfile', type=str, required=False, default=None, help='pre-existing file of variant locations in assembly compared to benchmark')
    parser.add_argument('--structureonly', action='store_true', required=False, help='analyse only the long-range structure of the assembly')
    parser.add_argument('-A', '--assembly', type=str, required=False, default="test", help='name of the assembly being tested--should be query sequence in bam file')
    parser.add_argument('-B', '--benchmark', type=str, required=False, default="truth", help='name of the assembly being used as a benchmark--should be the reference sequence in the bam file')
    parser.add_argument('-c', '--config', type=str, required=False, default="benchconfig.txt", help='path to a config file specifying locations of benchmark data files')
    parser.add_argument('--debug', action='store_true', required=False, help='print verbose output to log file for debugging purposes')

    return parser

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    if not args.bam and not args.paf:
        logger.critical("Must specify either a bam file with --bam or a paf file with --paf")
        exit(1)

    return args

def read_config_data(args)->dict:
    configfile = args.config
    configpath = Path(configfile)

    configvals = {}
    if not configpath.exists():
        template_res = importlib.resources.files("q100bench").joinpath('benchconfig.txt')
        logger.info("Using resource locations from default config file " + str(template_res))
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
        logger.info("Using resource locations from " + configfile)
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

    logfile = args.prefix + ".log"
    logformat = '%(asctime)s %(message)s'
    if args.debug:
        logging.basicConfig(filename=logfile, level=logging.DEBUG, format=logformat)
        logger.info('Logging verbose output for debugging.')
    else:
        logging.basicConfig(filename=logfile, level=logging.INFO, format=logformat)

    # check for necessary installed programs and write an output directory:
    check_for_bedtools()
    no_rscript = check_for_R()

    # dictionary of parameters from the benchmark configuration file:
    benchparams = read_config_data(args)

    # pysam objects for the benchmark and test assembly fasta files:
    ref = Path(args.reffasta)
    query = Path(args.queryfasta)
    if not ref.is_file() or not query.is_file():
        logger.critical("Ref fasta file " + args.reffasta + " and query fasta file " + args.queryfasta + " must exist and be readable")
        exit(1)
    refobj = pysam.FastaFile(args.reffasta)
    queryobj = pysam.FastaFile(args.queryfasta)

    outputdir = output.create_output_directory(args.prefix)

    # dictionary of this run's output file names:
    outputfiles = output.name_output_files(args, outputdir)

    logger.info("Step 1 (of 9): Writing bed files for excluded regions, test assembly scaffold spans, lengths, N stretches, and contigs (ignoring stretches of less than " + str(args.minns) + " Ns)")
    bedregiondict = {}
    # merged excluded regions are saved as BedTool object "allexcludedregions" in bedregiondict here:
    seqparse.write_genome_bedfiles(queryobj, refobj, args, benchparams, outputfiles, bedregiondict)
   
    # read in alignments from BAM or PAF format, filtering out secondaries and finding the longest increasing sequence along the reference:
    alignobj = None
    pafaligns = None
    if args.bam:
        if not os.access(args.bam, os.R_OK):
            logger.critical("BAM file " + args.bam + " is not readable")
            exit(1)
        alignobj = pysam.AlignmentFile(args.bam, "rb")
        aligndata = alignparse.read_bam_aligns(alignobj, args.minalignlength)
        rlis_aligndata = mummermethods.filter_aligns(aligndata, "target")
    else:
        pafaligns = alignparse.read_paf_aligns(args.paf, args.minalignlength)
        aligndata = pafaligns
        rlis_aligndata = mummermethods.filter_aligns(pafaligns, "target")

    # find general stats about contig/scaffold lengths, N/L50's, etc.:
    logger.info("Step 2 (of 9): Writing general statistics about " + args.assembly + " assembly")
    benchmark_stats = stats.write_general_assembly_stats(refobj, queryobj, bedregiondict["testnonnregions"], bedregiondict["testnregions"], outputfiles, args)

    # find clusters of consistent, covering alignments and calculate continuity statistics:
    logger.info("Step 3 (of 9): Assessing overall structural alignment of assembly")
    alignparse.assess_overall_structure(rlis_aligndata, refobj, queryobj, outputfiles, bedregiondict, benchmark_stats, args)
    structvar.write_structural_errors(refobj, queryobj, outputfiles, benchmark_stats, args)
    stats.write_aligned_cluster_stats(outputfiles, benchmark_stats, args)

    if not args.structureonly:
       logger.info("Step 4 (of 9): Writing bed files of regions covered by alignments of " + args.assembly + " to " + args.benchmark)
       # read in locations of het variants in the benchmark:
       hetsites = phasing.read_hetsites(benchparams["hetsitevariants"])
       hetarrays = phasing.sort_chrom_hetsite_arrays(hetsites)
       
       [refcoveredbed, querycoveredbed, variants, hetsitealleles] = alignparse.write_bedfiles(alignobj, pafaligns, refobj, queryobj, hetarrays, outputfiles["testmatcovered"], outputfiles["testpatcovered"], outputfiles["truthcovered"], outputfiles["coveredhetsitealleles"], bedregiondict["allexcludedregions"], args)

       # create merged unique outputfiles:
       [mergedtruthcoveredbed, outputfiles["mergedtruthcovered"]] = bedtoolslib.mergebed(outputfiles["truthcovered"])
       [mergedtestmatcoveredbed, outputfiles["mergedtestmatcovered"]] = bedtoolslib.mergebed(outputfiles["testmatcovered"])
       [mergedtestpatcoveredbed, outputfiles["mergedtestpatcovered"]] = bedtoolslib.mergebed(outputfiles["testpatcovered"])
   
       logger.info("Step 5 (of 9): Writing primary alignment statistics about " + args.assembly + " assembly")
       stats.write_merged_aligned_stats(refobj, queryobj, mergedtruthcoveredbed, mergedtestmatcoveredbed, mergedtestpatcoveredbed, outputfiles, benchmark_stats, args)

       if alignobj is not None:
   
           # classify variant errors as phasing or novel errors:
           logger.info("Step 6 (of 9): Writing phase switch statistics")
           stats.write_het_stats(outputfiles, benchmark_stats, args)
           logger.info("Step 7 (of 9): Determining whether errors are switched haplotype or novel")
           errors.classify_errors(refobj, queryobj, variants, hetsites, outputfiles, benchparams, benchmark_stats, args)
           stats.write_qv_stats(benchmark_stats, outputfiles, args)

           # evaluate mononucleotide runs:
           logger.info("Step 8 (of 9): Assessing accuracy of mononucleotide runs")
           bedtoolslib.intersectbed(benchparams["mononucruns"], outputfiles["mergedtruthcovered"], outputfile=outputfiles["coveredmononucsfile"], writefirst=True)
           mononucswithvariantsbedfile = bedtoolslib.intersectbed(outputfiles["coveredmononucsfile"], outputfiles["bencherrortypebed"], outputfiles["mononucswithvariantsfile"], outerjoin=True, writeboth=True)
           mononucstats = errors.gather_mononuc_stats(outputfiles["mononucswithvariantsfile"], outputfiles["mononucstatsfile"])
           stats.write_mononuc_stats(mononucstats, outputfiles, benchmark_stats, args)
    
    # plot alignment coverage across assembly and genome:
    if not no_rscript:
        logger.info("Step 9 (of 9): Creating plots")
        if not args.structureonly:
            plots.plot_benchmark_align_coverage(args.assembly, args.benchmark, outputdir, benchparams["resourcedir"])
            plots.plot_testassembly_align_coverage(args.assembly, outputdir, benchparams["resourcedir"])
            plots.plot_assembly_error_stats(args.assembly, args.benchmark, outputdir)
            if alignobj is not None:
                plots.plot_mononuc_accuracy(args.assembly, outputdir, benchparams["resourcedir"])
        plots.plot_svcluster_align_plots(args.assembly, args.benchmark, outputfiles["alignplotdir"], benchparams["resourcedir"], refobj)


if __name__ == "__main__":
    main()
