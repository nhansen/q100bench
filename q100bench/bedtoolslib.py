import os
import re
import pybedtools

def mergebed(bedfile:str)->str:

    mergedbed = bedfile.replace(".bed", ".merged.bed")

    unmergedints = pybedtools.bedtool.BedTool(bedfile)
    numunmerged = unmergedints.count()
   
    if numunmerged > 0:
        mergedints = unmergedints.merge(c=4, o='collapse', delim="|")
        mergedints.saveas(mergedbed)
    else:
        with open(mergedbed, 'a'):
            os.utime(mergedbed, None) 
        mergedints = unmergedints

    return [mergedints, mergedbed]

def intersectbed(bedfile1:str, bedfile2:str, outputfile:str, writefirst=False, writeboth=False, outerjoin=False)->str:

    print("Intersecting " + bedfile1 + " and " + bedfile2 + " to create " + outputfile)

    command = "bedtools intersect -a " + bedfile1 + " -b " + bedfile2
    if writefirst:
        command = command + " -wa"
    if writeboth:
        command = command + " -wo"
    if outerjoin:
        command = command + " -loj"
    os.system(command + " > " + outputfile)

    return outputfile

