import re
import sys
import pysam
import pybedtools
import logging
from pathlib import Path
from collections import namedtuple
from q100bench import seqparse

# create namedtuple for bed intervals:
varianttuple = namedtuple('varianttuple', ['chrom', 'start', 'end', 'name', 'vartype', 'excluded']) 

logger = logging.getLogger(__name__)

def read_hetsites(hetsitefile)->dict:
    hets = Path(hetsitefile)
    hetsites = {}
    if not hets.is_file():
        logger.error("Het site file " + hetsitefile + " does not exist so phasing analysis will be incomplete")
    else:
        refobj = pysam.FastaFile(args.reffasta)
        with open(hetsitefile, "r") as hfh:
            hetsiteline = hfh.readline()
            while hetsiteline:
                hetsiteline = hetsiteline.rstrip()
                [chrom, start, end, name] = hetsiteline.split("\t")
                namefields = name.split("_")
                refallele = namefields[-3]
                altallele = namefields[-2]
                hetsitename = chrom + "_" + str(int(start) + 1) + "_" + refallele + "_" + altallele 
                if refallele != "*" and altallele != "*" and len(refallele)==1 and len(altallele)==1:
                    vartype = 'SNV'
                else:
                    vartype = 'INDEL'
                hetsites[hetsitename] = varianttuple(chrom=chrom, start=int(start), end=int(end), name=name, vartype=vartype, excluded=False)
                hetsiteline = hfh.readline()

    return hetsites

def sort_chrom_hetsite_arrays(hetsites:dict):

    chromhetsites = {}
  
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
    logger.debug("Opening " + hetbed + " to write het alleles along assembly contigs")
    with open(hetbed, "w") as hfh:
        for contig in sorted(contigsortedhetalleles.keys()):
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
            


