import re
import sys
import random
import pybedtools
import logging
from collections import namedtuple
from q100bench import seqparse
from q100bench import alignparse
from q100bench import bedtoolslib

# create namedtuple for bed intervals:
varianttuple = namedtuple('varianttuple', ['chrom', 'start', 'end', 'name', 'vartype', 'excluded']) 

logger = logging.getLogger(__name__)

def classify_errors(refobj, queryobj, variants, hetsites, outputdict, benchparams, stats, args)->str:

    bencherrorfile = outputdict["bencherrortypebed"]
    stats["singlebasecounts"] = {}
    stats["indellengthcounts"] = {}
    stats["totalerrorsinaligns"] = 0

    testerrorfile = outputdict["testerrortypebed"]
    tfh = open(testerrorfile, "w")

    excludederrorfile = outputdict["benchexcludederrortypebed"]
    xfh = open(excludederrorfile, "w")

    with open(bencherrorfile, "w") as efh:
        for variant in variants:
            namefields = variant.name.split("_")
            numfields = len(namefields)
            contigname = "_".join(namefields[0:numfields-4])
            pos = int(namefields[-4])
            refallele = namefields[-3]
            altallele = namefields[-2]
            if namefields[-1] == "F":
                alignstrand = '+'
            else:
                alignstrand = '-'
            logger.debug("Splitting variant name " + variant.name + " and found contigname " + contigname)
            varname = variant.chrom + "_" + str(int(variant.start) + 1) + "_" + refallele + "_" + altallele

            if varname in hetsites.keys():
                errortype = 'PHASING'
                errortypecolor = '255,0,0'
            else:
                errortype = 'CONSENSUS'
                errortypecolor = '0,0,255'

            if variant.excluded:
                xfh.write(variant.chrom + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + varname + "\t1000\t" + alignstrand + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + errortypecolor + "\t" + errortype + "\t" + variant.vartype + "\t" + variant.name + "\n")
            else:
                efh.write(variant.chrom + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + varname + "\t1000\t" + alignstrand + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + errortypecolor + "\t" + errortype + "\t" + variant.vartype + "\t" + variant.name + "\n")

            tfh.write(contigname + "\t" + str(pos-1) + "\t" + str(pos - 1 + len(altallele)) + "\t" + variant.name + "\t1000\t" + alignstrand + "\t" + str(pos-1) + "\t" + str(pos - 1 + len(altallele)) + "\t" + errortypecolor + "\t" + errortype + "\t" + variant.vartype + "\t" + varname + "\n")
            stats["totalerrorsinaligns"] = stats["totalerrorsinaligns"] + 1
                
            # tally statistics for non-phasing errors:

            if variant.excluded or errortype != "CONSENSUS":
                continue
            if variant.vartype == "SNV":
                snvkey = refallele + "_" + altallele
                if snvkey in stats["singlebasecounts"]:
                    stats["singlebasecounts"][snvkey] = stats["singlebasecounts"][snvkey] + 1
                else:
                    stats["singlebasecounts"][snvkey] = 1
            else:
                trueref = refallele.replace('*', '')
                truealt = altallele.replace('*', '')
                lengthdiff = len(truealt) - len(trueref)
                if lengthdiff in stats["indellengthcounts"]:
                    stats["indellengthcounts"][lengthdiff] = stats["indellengthcounts"][lengthdiff] + 1
                else:
                    stats["indellengthcounts"][lengthdiff] = 1
    tfh.close()

    return stats

def gather_mononuc_stats(coveredmononucbedfile:str, mononucstatsfile:str):

    p = {}
    p["A"] = re.compile("^[aA]+$")
    p["T"] = re.compile("^[tT]+$")
    p["C"] = re.compile("^[cC]+$")
    p["G"] = re.compile("^[gG]+$")
    result = {}
    sfh = open(mononucstatsfile, "w")
    with open(coveredmononucbedfile, "r") as mfh:
        mononucline = mfh.readline()
        while mononucline:
            mononucline = mononucline.rstrip()
            [chrom, start, end, name, score, strand, widestart, wideend, color, chrom_b, start_b, end_b, variant_name, score_b, strand_b, widestart_b, wideend_b, color_b, error_type, variant_type, queryvariantname] = mononucline.split("\t")
            runlength = int(end) - int(start)
            namefields = name.split("_")
            repeatedbase = namefields[-1]
            if chrom_b == '.':
                result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':runlength, 'type':'CORRECT'}
                sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(runlength) + "\tCORRECT\n")
            else:
                variantnamefields = variant_name.split("_")
                refbases = variantnamefields[-2]
                altbases = variantnamefields[-1]
                if refbases == "*" and p[repeatedbase].match(altbases): # increased length
                    newlength = runlength + len(altbases)
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':error_type}
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + error_type + "\n")
                elif altbases == "*" and p[repeatedbase].match(refbases): # decreased length
                    newlength = runlength - len(refbases)
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':error_type}
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + error_type + "\n")
                elif p[repeatedbase].match(refbases) and p[repeatedbase].match(altbases): # expanded notation
                    newlength = len(altbases)
                    reflength = len(refbases)
                    if reflength != runlength:
                        logger.warning("Mononuc var reflength doesn\'t match benchmark run length at " + chrom + ":" + str(start) + "-" + str(end))
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + error_type + "\n")
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':error_type}
                else: # complex error
                    newlength = -1
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':error_type}
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + error_type + "\n")

            mononucline = mfh.readline()

    return result

def assess_mononuc_read_coverage(align_obj, mononucbedfile, outputdict, bedintervals, hetsitedict, args):

    p = {}
    p["A"] = re.compile("^[aA]+$")
    p["T"] = re.compile("^[tT]+$")
    p["C"] = re.compile("^[cC]+$")
    p["G"] = re.compile("^[gG]+$")

    mononucstatsfile = outputdict["mononucstatsfile"]

    mononucbedintervals = pybedtools.BedTool(mononucbedfile)
    if bedintervals is not None:
        includedmononucbeds = bedtoolslib.intersectintervals(mononucbedintervals, bedintervals, wa=True)
        includedmononucbeds.saveas(outputdict["includedmononucfile"])
    else:
        mononucbedintervals.saveas(outputdict["includedmononucfile"])

    mononucdict = {}
    msfh = open(mononucstatsfile, "w", buffering=1)
    # for each mononuc run in the benchmark, retrieve reads, count length of mononuc, classify as CORRECT, HET, ERROR, or COMPLEX
    with open(outputdict["includedmononucfile"], "r") as mfh:
        mononucline = mfh.readline()
        while mononucline:
            mononucline = mononucline.rstrip()
            [chrom, start, end, name, score, strand, widestart, wideend, color] = mononucline.split("\t")
            runlength = int(end) - int(start)
            if runlength not in mononucdict:
                mononucdict[runlength] = {}
            namefields = name.split("_")
            repeatedbase = namefields[-1]
            refalleleseq = repeatedbase * runlength
            for readalign in align_obj.fetch(contig=chrom, start=int(start), stop=int(end)):
                if args.downsample is not None and random.random() >= args.downsample:
                    continue
                queryseq = readalign.query_sequence
                if readalign.is_secondary or queryseq is None:
                    continue
                readname = readalign.query_name
                if readalign.reference_start >= int(start) or readalign.reference_end <= int(end):
                    continue
                pairs = readalign.get_aligned_pairs()
                if len(pairs) == 0:
                    continue
                # find zero-based read pos of base aligned to ref base one before mononuc:
                readstart = find_readpos_in_pairs(pairs, int(start)-1)
                # find zero-based read pos of base aligned to ref base one after mononuc:
                readend = find_readpos_in_pairs(pairs, int(end))
                if readstart is not None and readend is not None:
                    queryseq = queryseq.upper()
                    alleleseq = queryseq[readstart+1:readend]
                    if p[repeatedbase].match(alleleseq):
                        numbases = len(alleleseq)
                        matchtype = "CORRECT"
                        if numbases != runlength:
                            potentialhetname = chrom + "_" + str(int(start)+1) + "_" + refalleleseq + "_" + alleleseq
                            if potentialhetname in hetsitedict:
                                matchtype = "HET"
                            else:
                                matchtype = "ERROR"
                    else:
                        numbases = -1
                        matchtype = "COMPLEX"

                    msfh.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + readname + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(numbases) + "\t" + matchtype + "\n")
                    if numbases not in mononucdict[runlength]:
                        mononucdict[runlength][numbases] = {}
                    if matchtype not in mononucdict[runlength][numbases]:
                        mononucdict[runlength][numbases][matchtype] = 1
                    else:
                        mononucdict[runlength][numbases][matchtype] = mononucdict[runlength][numbases][matchtype] + 1
                else:
                    logger.debug(readname + " is unaligned at one endpoint!")
                    if readstart is None:
                        logger.warning(readname + " " + str(readalign.reference_start) + "-" + str(readalign.reference_end) + " does not have a start")
                    if readend is None:
                        logger.warning(readname + " " + str(readalign.reference_start) + "-" + str(readalign.reference_end) + " does not have a end")

            mononucline = mfh.readline()

    return mononucdict

# "pairs" is a structure created by pysam for an alignment, and contains a list
#   of tuples "consisting of the 0-based offset from the start of the read 
#   sequence followed by the 0-based reference position
#
# This routine searches the list of tuples for the one corresponding to the 
#   desired reference position ("pos"), and returns the read position
#   for that tuple, assuming it exists
#
# Why would a tuple have a ref position of "None"? These tuples are insertions
#   in the read sequence (and deletions from the reference have read pos of 
#   None). In the case of deletions from the reference, the "ifnone" argument
#   can tell the routine whether to report the read position that is one to 
#   the left of the deleted ref location ("lower") or one to the right ("higher")
#
# Note: this routine takes a *zero-based* position as its "pos" argument!
#
def find_readpos_in_pairs(pairs, pos, ifnone="lower"):
    # number of tuples:
    alignlength = len(pairs)
    ilow = 0
    ihigh = alignlength - 1

    hiref = pairs[ihigh][1]
    while hiref is None and ihigh > 0:
        ihigh = ihigh - 1
        hiref = pairs[ihigh][1]
    lowref = pairs[ilow][1]
    while lowref is None and ilow < len(pairs):
        ilow = ilow + 1
        lowref = pairs[ilow][1]

    # these shouldn't really happen
    if lowref is None or hiref is None or lowref > pos or hiref < pos or lowref > hiref:
        return None

    # need to be careful here to avoid an infinite loop:
    lastimid = None
    while ihigh - 1 > ilow and hiref >= pos and lowref <= pos:
        imid = int((ilow + ihigh)/2)
        midref = pairs[imid][1]
        while midref is None:
            if ifnone == "lower":
                imid = imid - 1
            else:
                imid = imid + 1
            midref = pairs[imid][1]
            if (ifnone == "lower" and imid <= ilow) or (ifnone=="higher" and imid >= ihigh):
                break
        if midref is None:
            return None
        if lastimid is not None and lastimid == imid:
            for i in range(ilow, ihigh):
                if pairs[i][1] is not None and pairs[i][1] == pos:
                    return pairs[i][0]
            return None
        
        if midref == pos:
            return pairs[imid][0] # might be None if nothing is aligned here
        elif midref > pos:
            ihigh = imid
            hiref = pairs[ihigh][1]
        elif midref < pos:
            ilow = imid
            lowref = pairs[ilow][1]
        if (pairs[ilow][1] is not None and pairs[ilow][1] > pos) or (pairs[ihigh][1] is not None and pairs[ihigh][1] < pos):
            return None
        lastimid = imid

    return None

def assess_read_align_errors(align_obj, refobj, readerrorfile:str, bedintervals, hetsitedict, args):
   
    if not args.errorfile:
        refh = open(readerrorfile, "w")

    stats = {}
    stats["totalalignedbases"] = 0
    stats["totalclippedbases"] = 0
    stats["totalerrorsinaligns"] = 0
    stats["singlebasecounts"] = {}
    stats["indellengthcounts"] = {}
    alignsprocessed = 0

    # are we dealing with just a set of regions? or the whole genome?
    if bedintervals is None:
        bedintervals = [None]

    for benchinterval in bedintervals:
        if benchinterval is not None:
            regionstring = benchinterval.chrom + ":" + str(benchinterval.start + 1) + "-" + str(benchinterval.end)
        else:
            regionstring = None

        for align in align_obj.fetch(region=regionstring):
            if align.is_secondary or align.cigartuples is None or (args.downsample is not None and random.random() >= args.downsample):
                continue
   
            stats["totalalignedbases"] = stats["totalalignedbases"] + align.reference_length
            if benchinterval is not None: # shorten aligned length contribution if necessary
                if benchinterval.start > align.reference_start:
                    stats["totalalignedbases"] = stats["totalalignedbases"] + align.reference_start - benchinterval.start
                if benchinterval.end < align.reference_end:
                    stats["totalalignedbases"] = stats["totalalignedbases"] - align.reference_end + benchinterval.end

            # TODO: clipping calc should be amended to not count if alignment endpoints are outside desired benchinterval!
            cigartuples = align.cigartuples
            if cigartuples[0][0] in [4, 5]:
                if benchinterval is None or benchinterval.start < align.reference_start:
                    stats["totalclippedbases"] = stats["totalclippedbases"] + cigartuples[0][1]
            if cigartuples[-1][0] in [4, 5]:
                if benchinterval is None or benchinterval.end > align.reference_end:
                    stats["totalclippedbases"] = stats["totalclippedbases"] + cigartuples[-1][1]
            
            if not args.errorfile:
                query, querystart, queryend, ref, refstart, refend, strand = alignparse.retrieve_align_data(align)
                if strand == "F":
                    queryleft = querystart
                    queryright = queryend
                else:
                    queryleft = queryend
                    queryright = querystart
    
                hetsites = {}
                hetsitealleles = {} # no need to track het site alleles in this context
                queryobj = None
    
                read_variants = alignparse.align_variants(align, queryobj, query, querystart, queryend, refobj, ref, refstart, refend, strand, hetsites, hetsitealleles, True)
    
                for variant in read_variants:
                    namefields = variant.name.split("_")
                    pos = int(namefields[-4])
                    refallele = namefields[-3]
                    altallele = namefields[-2]
                    posend = pos + len(refallele) # this will be one base after pos for insertions, one base past end of ref allele otherwise
                    if benchinterval is not None and (variant.start < benchinterval.start or variant.end > benchinterval.end):
                        continue
                    if namefields[-1] == "F":
                        alignstrand = '+'
                    else:
                        alignstrand = '-'
                    varname = variant.chrom + "_" + str(int(variant.start) + 1) + "_" + refallele + "_" + altallele
        
                    if varname in hetsitedict.keys():
                        errortype = 'HET'
                        errortypecolor = '255,0,0'
                    else:
                        errortype = 'ERROR'
                        errortypecolor = '0,0,255'
    
                    # record error in stats:
                    if errortype == 'ERROR':
                        stats["totalerrorsinaligns"] = stats["totalerrorsinaligns"] + 1
                        if refallele != "*" and altallele != "*" and len(refallele) == 1 and len(altallele) == 1:
                            if alignstrand == "+":
                                readrefallele = refallele
                                readaltallele = altallele
                            else:
                                readrefallele = seqparse.revcomp(refallele)
                                readaltallele = seqparse.revcomp(altallele)

                            snvkey = readrefallele + "_" + readaltallele
                            if snvkey in stats["singlebasecounts"]:
                                stats["singlebasecounts"][snvkey] = stats["singlebasecounts"][snvkey] + 1
                            else:
                                stats["singlebasecounts"][snvkey] = 1
                        else:
                            trueref = refallele.replace('*', '')
                            truealt = altallele.replace('*', '')
                            lengthdiff = len(truealt) - len(trueref)
                            if lengthdiff in stats["indellengthcounts"]:
                                stats["indellengthcounts"][lengthdiff] = stats["indellengthcounts"][lengthdiff] + 1
                            else:
                                stats["indellengthcounts"][lengthdiff] = 1
    
                    refh.write(variant.chrom + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + varname + "\t1000\t" + alignstrand + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + errortypecolor + "\t" + errortype + "\t" + variant.name + "\n")
    
            alignsprocessed = alignsprocessed + 1
            if alignsprocessed == 100000*int(alignsprocessed/100000):
                logger.debug("Processed " + str(alignsprocessed) + " aligns")

    if args.errorfile:
        with open(args.errorfile, "r") as efh:
            errorline = efh.readline()
            while errorline:
                errorline = errorline.rstrip()
                [chrom, start, end, varname, score, strand, widestart, wideend, color, variant_type, queryvariantname] = errorline.split("\t")
                if variant_type == "HET":
                    errorline = efh.readline()
                    continue

                namefields = varname.split("_")
                refallele = namefields[-2]
                altallele = namefields[-1]
                
                # record error in stats:
                stats["totalerrorsinaligns"] = stats["totalerrorsinaligns"] + 1
                if refallele != "*" and altallele != "*" and len(refallele) == 1 and len(altallele) == 1:
                    snvkey = refallele + "_" + altallele
                    if snvkey in stats["singlebasecounts"]:
                        stats["singlebasecounts"][snvkey] = stats["singlebasecounts"][snvkey] + 1
                    else:
                        stats["singlebasecounts"][snvkey] = 1
                else:
                    trueref = refallele.replace('*', '')
                    truealt = altallele.replace('*', '')
                    lengthdiff = len(truealt) - len(trueref)
                    if lengthdiff in stats["indellengthcounts"]:
                        stats["indellengthcounts"][lengthdiff] = stats["indellengthcounts"][lengthdiff] + 1
                    else:
                        stats["indellengthcounts"][lengthdiff] = 1
                errorline = efh.readline()

    return stats
