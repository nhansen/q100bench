import re
import sys
import random
from collections import namedtuple
from q100bench import alignparse

# create namedtuple for bed intervals:
bedinterval = namedtuple('bedinterval', ['chrom', 'start', 'end', 'name', 'rest']) 

def classify_errors(refobj, queryobj, variants, hetsites, outputdict, benchparams, stats, args)->str:

    bencherrorfile = outputdict["bencherrortypebed"]
    stats["singlebasecounts"] = {}
    stats["indellengthcounts"] = {}
    stats["totalerrorsinaligns"] = 0

    testerrorfile = outputdict["testerrortypebed"]
    tfh = open(testerrorfile, "w")

    with open(bencherrorfile, "w") as efh:
        for variant in variants:
            namefields = variant.name.split("_")
            numfields = len(namefields)
            contigname = "_".join(namefields[0:numfields-5])
            pos = int(namefields[-4])
            refallele = namefields[-3]
            altallele = namefields[-2]
            if namefields[-1] == "F":
                alignstrand = '+'
            else:
                alignstrand = '-'
            varname = variant.chrom + "_" + str(int(variant.start) + 1) + "_" + refallele + "_" + altallele

            if varname in hetsites.keys():
                errortype = 'PHASING'
                errortypecolor = '255,0,0'
            else:
                errortype = 'CONSENSUS'
                errortypecolor = '0,0,255'
    
            efh.write(variant.chrom + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + varname + "\t1000\t" + alignstrand + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + errortypecolor + "\t" + errortype + "\t" + variant.name + "\n")
            tfh.write(contigname + "\t" + str(pos-1) + "\t" + str(pos - 1 + len(altallele)) + "\t" + variant.name + "\t1000\t" + alignstrand + "\t" + str(pos-1) + "\t" + str(pos - 1 + len(altallele)) + "\t" + errortypecolor + "\t" + errortype + "\t" + varname + "\n")
            stats["totalerrorsinaligns"] = stats["totalerrorsinaligns"] + 1
                
            # tally statistics for non-phasing errors:

            if errortype != "CONSENSUS":
                continue
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
            [chrom, start, end, name, score, strand, widestart, wideend, color, chrom_b, start_b, end_b, variant_name, score_b, strand_b, widestart_b, wideend_b, color_b, variant_type, queryvariantname] = mononucline.split("\t")
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
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':variant_type}
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + variant_type + "\n")
                elif altbases == "*" and p[repeatedbase].match(refbases): # decreased length
                    newlength = runlength - len(refbases)
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':variant_type}
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + variant_type + "\n")
                elif p[repeatedbase].match(refbases) and p[repeatedbase].match(altbases): # expanded notation
                    newlength = len(altbases)
                    reflength = len(refbases)
                    if reflength != runlength:
                        print("Mononuc var reflength doesn\'t match benchmark run length at " + chrom + ":" + str(start) + "-" + str(end), file=sys.stderr)
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + variant_type + "\n")
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':variant_type}
                else: # complex error
                    newlength = -1
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':variant_type}
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + variant_type + "\n")

            mononucline = mfh.readline()

    return result

def assess_mononuc_read_coverage(align_obj, region:str, mononucbedfile:str, mononucstatsfile:str, hetsitedict, args):

    p = {}
    p["A"] = re.compile("^[aA]+$")
    p["T"] = re.compile("^[tT]+$")
    p["C"] = re.compile("^[cC]+$")
    p["G"] = re.compile("^[gG]+$")
  
    mononucdict = {}
    msfh = open(mononucstatsfile, "w", buffering=1)
    with open(mononucbedfile, "r") as mfh:
        mononucline = mfh.readline()
        while mononucline:
            mononucline = mononucline.rstrip()
            [chrom, start, end, name, score, strand, widestart, wideend, color] = mononucline.split("\t")
            runlength = int(end) - int(start)
            namefields = name.split("_")
            repeatedbase = namefields[-1]
            refalleleseq = repeatedbase * runlength
            #print(chrom + ":" + str(start) + "-" + str(end) + " " + repeatedbase + "\t" + str(runlength))
            for readalign in align_obj.fetch(contig=chrom, start=int(start), stop=int(end)):
                if args.downsample is not None and random.random() >= args.downsample:
                    continue
                pairs = readalign.get_aligned_pairs()
                readname = readalign.query_name
                if len(pairs) == 0:
                    continue
                readstart = find_readpos_in_pairs(pairs, int(start)-1)
                readend = find_readpos_in_pairs(pairs, int(end)+1)
                if readstart is not None and readend is not None:
                    queryseq = readalign.query_sequence
                    if readalign.is_secondary or queryseq is None:
                        continue
                    queryseq = queryseq.upper()
                    alleleseq = queryseq[readstart:readend]
                    alleleseq = alleleseq[1:-1]
                    if p[repeatedbase].match(alleleseq):
                        numbases = len(alleleseq)
                        type = "CORRECT"
                        if numbases != runlength:
                            potentialhetname = chrom + "_" + str(int(start)+1) + "_" + refalleleseq + "_" + alleleseq
                            if potentialhetname in hetsitedict:
                                type = "HET"
                            else:
                                type = "ERROR"
                        msfh.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + readname + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(numbases) + "\t" + type + "\n")
                    else:
                        msfh.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + readname + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + "COMPLEX" + "\t" + type + "\n")
                else:
                    print(readname + " is unaligned at one endpoint!")

            mononucline = mfh.readline()

    return 0

def find_readpos_in_pairs(pairs, pos):
    alignlength = len(pairs)
    low = 0
    high = alignlength - 1

    hiref = pairs[high][1]
    while hiref is None and high > 0:
        high = high - 1
        hiref = pairs[high][1]
    lowref = pairs[low][1]
    while lowref is None and low < len(pairs):
        low = low + 1
        lowref = pairs[low][1]

    if lowref is None or hiref is None or lowref > pos or hiref < pos:
        return None

    while high > low and hiref >= pos and lowref <= pos:
        mid = int((low + high)/2)
        if random.random() >= 0.5:
            mid = mid + 1
        himid = mid
        lowmid = mid
        while himid < high and pairs[himid][1] is None:
            himid = himid + 1
        while lowmid > low and pairs[lowmid][1] is None:
            lowmid = lowmid - 1
        if pairs[lowmid][1] is None or pairs[himid][1] is None:
            return None
        if abs(pairs[lowmid][1] - pos) < abs(pairs[himid][1] - pos):
            mid = lowmid
        else:
            mid = himid
        midval = pairs[mid][1]
        #print(str(low) + " " + str(high) + " " + str(midval) + " " + str(pos) + " " + str(pairs[low][1]) + " " + str(pairs[high][1]) + " " + str(len(pairs)))
        if midval == pos:
            return pairs[mid][0] # might be None if nothing is aligned here
        elif midval > pos:
            high = mid
        elif midval < pos:
            low = mid
        if (pairs[low][1] is not None and pairs[low][1] > pos) or (pairs[high][1] is not None and pairs[high][1] < pos):
            return None

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
            if benchinterval is not None: # shorted aligned length contribution if necessary
                if benchinterval.start > align.reference_start:
                    #print("Shortening aligned read length for " + align.query_name + " by " + str(benchinterval.start - align.reference_start))
                    stats["totalalignedbases"] = stats["totalalignedbases"] + align.reference_start - benchinterval.start
                if benchinterval.end < align.reference_end:
                    #print("Shortening aligned read length for " + align.query_name + " by " + str(align.reference_end - benchinterval.end))
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
                query, querystart, queryend, ref, refstart, refend, strand = alignparse.retrieve_align_data(align, args)
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
                        #print("Skipping variant from " + str(variant.start) + " to " + str(variant.end) + " because not in interval " + regionstring)
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
    
                    refh.write(variant.chrom + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + varname + "\t1000\t" + alignstrand + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + errortypecolor + "\t" + errortype + "\t" + variant.name + "\n")
    
            alignsprocessed = alignsprocessed + 1
            if alignsprocessed == 100000*int(alignsprocessed/100000):
                print("Processed " + str(alignsprocessed) + " aligns")

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
