import re
import pysam
import pybedtools
import logging
from pathlib import Path
from q100bench import bedtoolslib
from collections import namedtuple
import time

logger = logging.getLogger(__name__)

def write_genome_bedfiles(queryobj, refobj, args, benchparams, outputfiles, bedobjects):

    write_excluded_bedfile(args, benchparams, outputfiles, bedobjects)
    write_test_genome_bedfile(queryobj, args, outputfiles, bedobjects)
    find_all_ns(queryobj, args, outputfiles, bedobjects)

def write_excluded_bedfile(args, benchparams, outputfiles, bedobjects):

    allexcludedbedfiles = []
    if args.excludefile is not None:
        allexcludedbedfiles.append(args.excludefile)
    if "excluderegions" in benchparams.keys():
        configexcluderegions = benchparams["excluderegions"]
        configexcludepath = Path(configexcluderegions)
        if configexcludepath.is_file():
            allexcludedbedfiles.append(configexcluderegions)
        else:
            logger.warning("Exclude file " + configexcluderegions + " does not exist so will not be used")
    if "nstretchregions" in benchparams.keys():
        nexcluderegions = benchparams["nstretchregions"]
        nexcludepath = Path(nexcluderegions)
        if nexcludepath.is_file():
            allexcludedbedfiles.append(nexcluderegions)
        else:
            logger.warning("Benchmark N file " + nexcluderegions + " does not exist so will not be used")
    if len(allexcludedbedfiles) == 0:
        bedobjects["allexcludedregions"] = pybedtools.BedTool("", from_string = True)
    elif len(allexcludedbedfiles) == 1:
        bedobjects["allexcludedregions"] = pybedtools.BedTool(allexcludedbedfiles[0])
    elif len(allexcludedbedfiles) > 1:
        logger.info("Merging " + str(allexcludedbedfiles) + " to create a combined excluded regions file")
        bedobjects["allexcludedregions"] = bedtoolslib.mergemultiplebedfiles(allexcludedbedfiles)
    bedobjects["allexcludedregions"].saveas(outputfiles["allexcludedbed"])

    return 0

def write_test_genome_bedfile(queryobj, args, outputfiles, bedobjects):

    genomebedstring = ""
    for scaffold in queryobj.references:
        scaffoldlength = queryobj.get_reference_length(scaffold)
        scaffstring = scaffold + "\t0\t" + str(scaffoldlength - 1) + "\n"
        genomebedstring += scaffstring
    bedobjects["testgenomeregions"] = pybedtools.BedTool(genomebedstring, from_string = True)
    bedobjects["testgenomeregions"].saveas(outputfiles["testgenomebed"])

    return 0

def find_all_ns(queryobj, args, outputfiles, bedobjects)->list:

    # did the command-line arguments specify a pre-existing bedfile of N-stretch locations?
    user_n_file = args.n_bedfile

    if not user_n_file:
        p = re.compile("N+")

    gapbedstring = ""
    contigbedstring = ""
    n_interval_dict = {}
    # if an n interval file has been specified in the command line options, read it into gapbedstring, and subtract from genome to get contigbedstring:
    if user_n_file:
        nbedobj = BedTool(user_n_file)
        for n_interval in nbedobj:
            chrom = n_interval.chrom
            start = n_interval.start
            end = n_interval.end
            gapname = n_interval.name
            gapbedstring += chrom + "\t" + str(start) + "\t" + str(end) + "\t" + gapname + "\n"
            if chrom in n_interval_dict.keys():
                n_interval_dict[chrom].append(n_interval)
            else:
                n_interval_dict[chrom] = [n_interval]
        for ref in queryobj.references:
            if ref not in n_interval_dict.keys():
                n_interval_dict[ref] = []
            contignum = 1
            contigstart = 0
            for interval in n_interval_dict[ref]:
                interval_name = ref + "." + str(contignum)
                contigbedstring += ref + "\t" + str(contigstart) + "\t" + str(interval.start) + "\t" + interval_name + "\n"
                contignum = contignum + 1
                contigstart = interval.end
            refend = queryobj.get_reference_length(ref)
            interval_name = ref + "." + str(contignum)
            contigbedstring += ref + "\t" + str(contigstart) + "\t" + str(refend) + "\t" + interval_name + "\n"
    else:
        for ref in queryobj.references:
            contignum = 1
            contigstart = 0
            findstring = 'N' * args.minns
            chromseq = queryobj.fetch(ref).upper()
            refend = queryobj.get_reference_length(ref)
            start = chromseq.find(findstring)
            while start != -1:
                end = start + args.minns
                while end < refend and chromseq[end] == 'N':
                    end = end + 1
                if end - start >= args.minns:
                    gapname = "N." + ref + "." + str(contignum)
                    gapstring = ref + "\t" + str(start) + "\t" + str(end) + "\t" + gapname + "\n"
                    gapbedstring += gapstring
                    contigend = start
                    contigname = ref + "." + str(contignum)
                    contigbedstring += ref + "\t" + str(contigstart) + "\t" + str(contigend) + "\t" + contigname + "\n"
                    contignum = contignum + 1
                    contigstart = end
                start = chromseq.find(findstring, end, refend)

            contigname = ref + "." + str(contignum)
            contigbedstring += ref + "\t" + str(contigstart) + "\t" + str(refend) + "\t" + contigname + "\n"

        bedobjects["testnonnregions"] = pybedtools.BedTool(contigbedstring, from_string = True)
        bedobjects["testnregions"] = pybedtools.BedTool(gapbedstring, from_string = True)

    if outputfiles["testnonnbed"]:
        bedobjects["testnonnregions"].saveas(outputfiles["testnonnbed"])
    if outputfiles["testnbed"]:
        bedobjects["testnregions"].saveas(outputfiles["testnbed"])

    return 0

def revcomp(seq:str) -> str:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'M': 'N', 'K': 'N', 'R': 'N', 'W': 'N', 'Y': 'N'}
    bases = list(seq)
    bases = bases[::-1]
    bases = [complement[base] for base in bases]

    return ''.join(bases)

