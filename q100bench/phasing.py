import re
import sys
import pysam
import pybedtools
import logging
from pathlib import Path
from collections import namedtuple
from q100bench import seqparse

# create namedtuple for bed intervals:
varianttuple = namedtuple('varianttuple', ['chrom', 'start', 'end', 'name', 'vartype', 'excluded', 'qvscore']) 

logger = logging.getLogger(__name__)

def read_hetsites(hetsitefile)->dict:
    hets = Path(hetsitefile)
    hetsites = {}
    if not hets.is_file():
        logger.error("Het site file " + hetsitefile + " does not exist so phasing analysis will be incomplete")
    else:
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
                hetsites[hetsitename] = varianttuple(chrom=chrom, start=int(start), end=int(end), name=name, vartype=vartype, excluded=False, qvscore=None)
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
    
# This routine duplicates the logic of "MerToPhaseBlock" from Arang Rhie's merqury package (but ported from Java to python)
# I've changed some variable names (e.g., "noBreak" to "includegaps") to make it easier (for me) to read:

def find_phase_blocks_from_marker_bed(bedfile:str, shortnum=1, shortlimit=200, includegaps=True)->dict:

    # determine two haplotypes in the "name" field of the input bed file (assumes only two distinct values plus optionally "gap")
    with open(bedfile, "r") as bfh:
        markerline = bfh.readline()
        hap1 = None
        hap2 = None
        while markerline:
            markerline = markerline.rstrip()
            [scaffold, start, end, markertype] = markerline.split("\t")
            if markertype != "gap" and hap1 is None:
                hap1 = markertype
            elif markertype != "gap" and markertype != hap1:
                hap2 = markertype
                break
            markerline = bfh.readline()
    if hap1 is not None and hap2 is not None:
        logger.info("Found " + hap1 + " and " + hap2 + " as haplotypes.")
        print("Found " + hap1 + " and " + hap2 + " as haplotypes.")
    else:
        logger.info("Unable to find two haplotypes in " + bedfile)
        print("Unable to find two haplotypes in " + bedfile)
        return {'phaseblocks':[], 'switches':[]}

    # find scaffold, endpoints, haplotype, and numbers of switches and markers for each phase block
    # also catalog switches (including "Same") in an array to return
    phaseblocks = []
    switches = []

    nummarkers = 0
    numswitches = 0
    numshortswitches = 0
    blockscaff = ""
    blockmarker = ""
    blockstart = -1
    blockend = -1
    shortstart = -1
    shortend = -1
    isshortswitch = True
    with open(bedfile, "r") as bfh:
        markerline = bfh.readline()
        while markerline:
            markerline = markerline.rstrip()
            [scaffold, start, end, markertype] = markerline.split("\t")
            start = int(start)
            end = int(end)

            if markertype == "gap": # gap
                if not includegaps:
                    # create a block with ending at the starting point of the gap, and add it to the list of phaseblocks:
                    blockend = start
                    if blockscaff != "":
                        print("Appending block " + blockscaff + ":" + str(blockstart) + "-" + str(blockend) + " type " + blockmarker + " with " + str(numswitches) + " switches and " + str(nummarkers) + " markers")
                        phaseblocks.append({'scaffold':blockscaff, 'blockstart':blockstart, 'blockend':blockend, 'blockmarker':blockmarker, 'numswitches':numswitches, 'nummarkers':nummarkers })
                    # start the next block at the end of the gap with "unknown" marker type, without switches and zero markers:
                    blockscaff = scaffold
                    blockstart = end
                    blockend = end
                    blockmarker = "unknown"
                    numswitches = 0
                    nummarkers = 0
                    isshortswitch = False
            elif (blockscaff == "" or blockscaff != scaffold): # scaffold has changed or this is start of the first scaff (not a gap line)
                # unless this is the start of the first scaffold, calculate number of markers in the block and add the current block to the list
                if blockscaff != "": # we are currently in a block
                    nummarkers = nummarkers - numshortswitches
                    print("Appending block " + blockscaff + ":" + str(blockstart) + "-" + str(blockend) + " type " + blockmarker + " with " + str(numswitches) + " switches and " + str(nummarkers) + " markers")
                    phaseblocks.append({'scaffold':blockscaff, 'blockstart':blockstart, 'blockend':blockend, 'blockmarker':blockmarker, 'numswitches':numswitches, 'nummarkers':nummarkers })
                    # short switches at the ends of scaffolds are included as long range switches
                    if numshortswitches > 0:
                        blockstart = shortstart
                        blockend = shortend
                        nummarkers = numshortswitches
                        blockmarker = markertype
                        print("Appending block " + blockscaff + ":" + str(blockstart) + "-" + str(blockend) + " type " + blockmarker + " with " + str(numswitches) + " switches and " + str(nummarkers) + " markers")
                        phaseblocks.append({'scaffold':blockscaff, 'blockstart':blockstart, 'blockend':blockend, 'blockmarker':blockmarker, 'numswitches':numswitches, 'nummarkers':nummarkers })
                # start the first block on this new scaffold
                blockmarker = markertype
                numswitches = 0
                numshortswitches = 0
                distfromswitch = 0
                nummarkers = 0
                blockscaff = scaffold
                blockstart = start
                blockend = end
                if (markertype == hap1 or markertype == hap2):
                    nummarkers = nummarkers + 1
                # java version writes out a switch with "Same" here
                print("SWITCH " + scaffold + ":" + str(start) + "-" + str(end) + " type " + markertype + " switch Same " + " previous type " + blockmarker)
                switches.append({'scaffold':scaffold, 'start':start, 'end':end, 'markertype':markertype, 'switchtype':'Same', 'blockmarker':blockmarker})
                isshortswitch = False
            else: # same scaffold, not a gap
                nummarkers = nummarkers + 1
                # if the marker type matches the block's type (or the last was an unknown at the beginning of a scaffold), extend block end:
                if blockmarker == markertype or blockmarker == "unknown":
                    blockend = end
                    if blockmarker == "unknown":
                        blockmarker = markertype
                    # if we're in a short switch, weve now come back to the block marker type, so reset shortswitch counter, but first add value to number of switches
                    if isshortswitch:
                        numswitches = numswitches + numshortswitches
                        numshortswitches = 0
                        isshortswitch = False
                        # java version writes out a switch with "SwitchBack" here
                        print("SWITCH " + scaffold + ":" + str(start) + "-" + str(end) + " type " + markertype + " switch SwitchBack " + " previous type " + blockmarker)
                        switches.append({'scaffold':scaffold, 'start':start, 'end':end, 'markertype':markertype, 'switchtype':'SwitchBack', 'blockmarker':blockmarker})
                    else:
                        # java version writes out a switch with "Same" here
                        print("SWITCH " + scaffold + ":" + str(start) + "-" + str(end) + " type " + markertype + " switch Same " + " previous type " + blockmarker)
                        switches.append({'scaffold':scaffold, 'start':start, 'end':end, 'markertype':markertype, 'switchtype':'Same', 'blockmarker':blockmarker})
                # if the marker is the opposite of the block's type, extend the short switch and increase the counter
                elif blockmarker != markertype:
                    shortend = end
                    numshortswitches = numshortswitches + 1
                    # if we were not in a short switch already, set the start of the switch region to this marker's start:
                    if not isshortswitch:
                        shortstart = start
                        isshortswitch = True
                        # java version writes out a switch with "Short" here
                        print("SWITCH " + scaffold + ":" + str(start) + "-" + str(end) + " type " + markertype + " switch Short " + " previous type " + blockmarker)
                        switches.append({'scaffold':scaffold, 'start':start, 'end':end, 'markertype':markertype, 'switchtype':'Short', 'blockmarker':blockmarker})
                    # if we were already in a short switch, check to see if it's numerous or long enough to be a long switch:
                    else:
                        distfromswitch = shortend - shortstart
                        # if we're in a long switch, decrease the number of markers by the number of short switches, and add a block
                        if numshortswitches >= shortnum or distfromswitch > shortlimit:
                            nummarkers = nummarkers - numshortswitches
                            print("Appending block " + scaffold + ":" + str(blockstart) + "-" + str(blockend) + " type " + blockmarker + " with " + str(numswitches) + " switches and " + str(nummarkers) + " markers")
                            phaseblocks.append({'scaffold':scaffold, 'blockstart':blockstart, 'blockend':blockend, 'blockmarker':blockmarker, 'numswitches':numswitches, 'nummarkers':nummarkers})
                            # then initiate a new block on the switched haplotype (allowing the formerly short block to extend)
                            blockstart = shortstart
                            blockend = end
                            nummarkers = numshortswitches
                            numswitches = 0
                            numshortswitches = 0
                            blockmarker = markertype
                            isshortswitch = False
                            # java version writes out a switch with "Long" here
                            print("SWITCH " + scaffold + ":" + str(start) + "-" + str(end) + " type " + markertype + " switch Long " + " previous type " + blockmarker)
                            switches.append({'scaffold':scaffold, 'start':start, 'end':end, 'markertype':markertype, 'switchtype':'Long', 'blockmarker':blockmarker})
            markerline = bfh.readline()

        print("Appending block " + scaffold + ":" + str(blockstart) + "-" + str(end) + " type " + blockmarker + " with " + str(numswitches) + " switches and " + str(nummarkers) + " markers")
        phaseblocks.append({'scaffold':scaffold, 'blockstart':blockstart, 'blockend':end, 'blockmarker':blockmarker, 'numswitches':numswitches,'nummarkers':nummarkers})

    return phaseblocks

