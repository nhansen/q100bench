import re
import sys
import pysam
import pybedtools
from collections import namedtuple
from q100bench import seqparse

# create namedtuple for bed intervals:
bedinterval = namedtuple('bedinterval', ['chrom', 'start', 'end', 'name', 'rest']) 

def read_hetsites(hetsitefile)->dict:
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
            hetsites[hetsitename] = bedinterval(chrom=chrom, start=int(start), end=int(end), name=name, rest='')
            hetsiteline = hfh.readline()

    return hetsites

def sort_chrom_hetsite_arrays(hetsites:dict):

    chromhetsites = {}
  
    #if hetsites.__class__==dict:
        #print("Hetsites is a dict!")
    #else:
        #print(str(hetsites.__class__), str(hetsites.__class__==dict))
    for hetsite in hetsites.values():
        if hetsite.__class__==dict: #dict
            hetsitetype = 'dict'
            chrom = hetsite['chrom']
            if chrom not in chromhetsites:
                chromhetsites[chrom] = []
            chromhetsites[chrom].append(hetsite)
        else: # tuple
            hetsitetype = 'tuple'
            chrom = hetsite.chrom
            if chrom not in chromhetsites:
                chromhetsites[chrom] = []
            chromhetsites[chrom].append(hetsite)
            
    for chrom in chromhetsites:
        if hetsitetype=="dict":
            chromhetsites[chrom].sort(key=lambda h: (h['start'], h['end']))
        else:
            chromhetsites[chrom].sort(key=lambda h: (h.start, h.end))
    
    return chromhetsites

def write_hetallele_bed(hetsitealleles:dict, hetbed:str):

    contigsortedhetalleles = sort_chrom_hetsite_arrays(hetsitealleles)
    print("Opening " + hetbed + " to write het alleles along assembly contigs")
    with open(hetbed, "w") as hfh:
        for contig in sorted(contigsortedhetalleles.keys()):
            numhets = len(contigsortedhetalleles)
            print("Processing " + str(numhets) + " hets for chrom " + contig)
            for hetsite in contigsortedhetalleles[contig]:
                hetname = hetsite['name']
                fields = hetname.split("_")
                strand = fields[-1]
                altallele = fields[-2]
                refallele = fields[-3]
                if hetsite['allele'] == refallele:
                    allelehap = 'SAMEHAP'
                elif hetsite['allele'] == altallele:
                    allelehap = 'ALTHAP'
                else:
                    allelehap = 'OTHER'
                assemblycontig = hetsite['query']
                assemblystart = hetsite['start'] - 1
                assemblyend = hetsite['end'] - 1
                hfh.write(contig + "\t" + str(hetsite['start']) + "\t" + str(hetsite['end']) + "\t" + hetsite['name'] + "\t" + hetsite['allele'] + "\t" + hetsite['ref'] + "\t" + str(hetsite['refstart']) + "\t" + str(hetsite['refend']) + "\t" + assemblycontig + "\t" + str(assemblystart) + "\t" + str(assemblyend) + "\t" + allelehap + "\n")
            


