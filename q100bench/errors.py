import re
import sys
from collections import namedtuple

# create namedtuple for bed intervals:
bedinterval = namedtuple('bedinterval', ['chrom', 'start', 'end', 'name', 'rest']) 

def classify_errors(refobj, queryobj, refcovered, querycovered, variants, outputdict, benchparams, args)->str:
    hetsitefile = benchparams["hetsitevariants"]

    hetsites = {}
    with open(hetsitefile, "r") as hfh:
        hetsiteline = hfh.readline()
        while hetsiteline:
            hetsiteline = hetsiteline.rstrip()
            [chrom, start, end, name] = hetsiteline.split("\t")
            namefields = name.split("_")
            refallele = namefields[-3]
            altallele = namefields[-2]
            hetsitename = chrom + "_" + str(int(start) + 1) + "_" + refallele + "_" + altallele 
            hetsites[hetsitename] = bedinterval(chrom=chrom, start=str(start), end=str(end), name=name, rest='')
            hetsiteline = hfh.readline()

    bencherrorfile = outputdict["bencherrortypebed"]

    testerrorfile = outputdict["testerrortypebed"]
    tfh = open(testerrorfile, "w")

    with open(bencherrorfile, "w") as efh:
        for variant in variants:
            namefields = variant.name.split("_")
            pos = int(namefields[-4])
            refallele = namefields[-3]
            altallele = namefields[-2]
            alignstrand = namefields[-1]
            varname = variant.chrom + "_" + str(int(variant.start) + 1) + "_" + refallele + "_" + altallele

            if varname in hetsites.keys():
                errortype = 'PHASING'
            else:
                errortype = 'CONSENSUS'
    
            efh.write(variant.chrom + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + varname + "\t" + errortype + "\n")
            tfh.write(str(pos-1) + "\t" + str(pos - 1 + len(altallele)) + "\t" + errortype + "\n")
    tfh.close()

    return bencherrorfile

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
            [chrom, start, end, name, score, strand, widestart, wideend, color, chrom_b, start_b, end_b, variant_name, variant_type] = mononucline.split("\t")
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

