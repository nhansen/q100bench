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

def mergeintervals(intervals):

    mergedints = intervals.merge(c=4, o='collapse', delim="|")

    return mergedints

def mergemultiplebedfiles(bedfilelist:list):

    if len(bedfilelist) < 2:
        print("Cannot call mergemultiplebedfiles on less than two bed files!")
        exit(1)

    firstbedtool = pybedtools.bedtool.BedTool(bedfilelist[0])
    allbedtools = firstbedtool.cat(bedfilelist[1])

    return allbedtools

def bedsum(intervals)->int:

    alllengths = map(len, intervals)
    bedsum = sum(alllengths)

    return bedsum

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
    intersectbed = pybedtools.bedtool.BedTool(outputfile)

    return [intersectbed, outputfile]

def intersectintervals(intervals1, intervals2, v=False, wa=False):

    intersectedints = intervals1.intersect(intervals2, v=v, wa=wa)

    return intersectedints

