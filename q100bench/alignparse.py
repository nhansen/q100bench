import re
import pysam
import pybedtools
from collections import namedtuple
from q100bench import seqparse

# create namedtuple for bed intervals:
bedinterval = namedtuple('bedinterval', ['chrom', 'start', 'end', 'name', 'rest']) 

# query start is always less than query end regardless of strand. query left always corresponds to ref start, 
# and so will be greater than query right when strand is reversed--all four are 1-based
def write_bedfiles(bamobj, refobj, queryobj, testmatbed, testpatbed, truthbed, variantbed, args):

    queryintervals = []
    refintervals = []
    user_variantfile = args.variantfile
    variants = []

    for align in bamobj.fetch():
        if align.is_secondary:
            continue
        if align.reference_length >= args.minalignlength:
            query, querystart, queryend, ref, refstart, refend, strand = retrieve_align_data(align, args)
            if strand == "F":
                queryleft = querystart
                queryright = queryend
            else:
                queryleft = queryend
                queryright = querystart
            querynamestring = query + "." + str(queryleft) + "." + str(queryright)
            refnamestring = ref + "." + str(refstart) + "." + str(refend) + "." + strand
            #print(querynamestring + ", " + refnamestring)
            queryintervals.append(bedinterval(chrom=query, start=querystart - 1, end=queryend, name=refnamestring, rest=''))
            refintervals.append(bedinterval(chrom=ref, start=refstart - 1, end=refend, name=querynamestring, rest=''))
            if user_variantfile is None:
                variants.extend(align_variants(align, queryobj, query, querystart, queryend, refobj, ref, refstart, refend, strand))

    phap1 = re.compile(r'.*MAT.*')
    phap2 = re.compile(r'.*PAT.*')
    with open(testmatbed, "w") as tmb:
        for testtuple in sorted(queryintervals, key=lambda h: (h.chrom, h.start, h.end)):
            if phap1.match(testtuple.name):
                tmb.write(testtuple.chrom + "\t" + str(testtuple.start) + "\t" + str(testtuple.end) + "\t" + testtuple.name + "\n")

    with open(testpatbed, "w") as tpb:
        for testtuple in sorted(queryintervals, key=lambda h: (h.chrom, h.start, h.end)):
            if phap2.match(testtuple.name):
                tpb.write(testtuple.chrom + "\t" + str(testtuple.start) + "\t" + str(testtuple.end) + "\t" + testtuple.name + "\n")

    with open(truthbed, "w") as rb:
        for truthtuple in sorted(refintervals, key=lambda h: (h.chrom, h.start, h.end)):
            rb.write(truthtuple.chrom + "\t" + str(truthtuple.start) + "\t" + str(truthtuple.end) + "\t" + truthtuple.name + "\n")

    if user_variantfile is None:
        with open(variantbed, "w") as vb:
            for varianttuple in sorted(variants, key=lambda h: (h.chrom, h.start, h.end)):
                vb.write(varianttuple.chrom + "\t" + str(varianttuple.start) + "\t" + str(varianttuple.end) + "\t" + varianttuple.name + "\n")

    if user_variantfile is not None:
        with open(user_variantfile, "r") as vh:
            variantline = vh.readline()
            while variantline:
                variantline = variantline.rstrip()
                [chrom, start, end, name] = variantline.split("\t")
                variants.append(bedinterval(chrom=chrom, start=str(start), end=str(end), name=name, rest=''))
                variantline = vh.readline()

    return [refintervals, queryintervals, variants]

def retrieve_align_data(align, args)->list:
    if align.is_reverse:
        strand = 'R'
    else:
        strand = 'F'
    query = align.query_name

    # number of hard clipped bases:
    cigartuples = align.cigartuples
    # will use actual one-based positions, so I don't go crazy, then report BED format zero-based half open
    if strand == 'F':
        if cigartuples[0][0] == 5:
            hardclip = cigartuples[0][1]
        else:
            hardclip = 0
        querystart = hardclip + align.query_alignment_start + 1 # first unclipped base (was 0-based, now 1-based)
        queryend = hardclip + align.query_alignment_end # was 0-based index just past the last unclipped base now 1-based last unclipped base
    elif strand == 'R':
        if cigartuples[-1][0] == 5:
            hardclip = cigartuples[-1][1]
        else:
            hardclip = 0
        querystart = hardclip + align.query_length - align.query_alignment_end + 1 # was 0-based index just past end of reversed query, now 1-based 5' end of query
        queryend = hardclip + align.query_length - align.query_alignment_start # was 0-based index of first unclipped base, now 1-based 3' unclipped end of query

    ref = align.reference_name

    refstart = align.reference_start + 1 # was zero-based, now one-based
    refend = align.reference_end #  was zero-based, adding one for bed format "half open", now one-based

    aligndata = [query, querystart, queryend, ref, refstart, refend, strand]

    return aligndata

# query start and query end are the lower and higher endpoints of the query seq in query coordinates (1-based)
# regardless of orientation of the alignment
def align_variants(align, queryobj, query:str, querystart:int, queryend:int, refobj, ref:str, refstart:int, refend:int, strand:str, widen=True)->list:

    #print("In align_variants with " + query + ":" + str(querystart) + "-" + str(queryend) + " " + ref + ":" + str(refstart) + "-" + str(refend) + "/" + strand)
    variantlist = []
    coveredregionlist = []
    homozygousregionlist = []

    queryseq = queryobj.fetch(reference=query, start=querystart-1, end=queryend).upper()
    refseq = refobj.fetch(reference=ref, start=refstart-1, end=refend).upper()

    alignstring = query + ":" + str(querystart-1) + "-" + str(queryend) + ";" + query + ":" + str(querystart-1) + "-" + str(queryend)
    strandsign = 1
    bedstrand = '+'
    if strand == 'R':
        queryseq = seqparse.revcomp(queryseq)
        strandsign = -1
        bedstrand = '-'

    alignops = align.cigartuples

    # first position will begin at 5'-most ref base of the alignment (regardless of strand)
    refcurrentoffset = 0
    querycurrentoffset = 0

    alignopindex = 0
    reflength = refend - refstart + 1
    querylength = queryend - querystart + 1
    while refcurrentoffset <= refend-refstart and alignopindex < len(alignops):
        alignop = alignops[alignopindex]
        op = alignop[0]
        oplength = alignop[1]

        if op in [0, 7, 8]: # MX= find SNV and MNVs
            for blockoffset in range(oplength):
                refpos = refcurrentoffset + blockoffset
                querypos = querycurrentoffset + blockoffset
                if refseq[refpos] != queryseq[querypos]:
                    if strand == 'F':
                        querycoordinate = querypos + querystart
                    else:
                        querycoordinate = queryend - querypos
                    variantname=query+"_"+str(querycoordinate)+"_"+refseq[refpos]+"_"+queryseq[querypos]+"_"+strand
                    additionalfields = "0\t" + bedstrand + "\t" + str(refpos+refstart-1) + "\t" + str(refpos+refstart) + "\t0,0,0\t" + alignstring
                    variantlist.append(bedinterval(chrom=ref, start=refpos+refstart-1, end=refpos+refstart, name=variantname, rest=additionalfields))

        if op in [2, 3]:
            refallele = refseq[refcurrentoffset:refcurrentoffset+oplength]
            queryallele = "*"
            if widen is True:
                extend = 0
                while querycurrentoffset + extend < querylength and refcurrentoffset + extend < reflength and refseq[refcurrentoffset + extend] == queryseq[querycurrentoffset + extend]:
                    #print(str(extend) + ": " + queryseq[querycurrentoffset + extend])
                    refallele = refallele + queryseq[querycurrentoffset + extend]
                    if queryallele == "*":
                        queryallele = queryseq[querycurrentoffset + extend]
                    else:
                        queryallele = queryallele + queryseq[querycurrentoffset + extend]
                    extend = extend + 1
            if strand == 'F':
                querycoordinate = querystart + querycurrentoffset
            else:
                querycoordinate = queryend - querycurrentoffset
            variantname=query+"_"+str(querycoordinate)+"_"+refallele+"_"+queryallele+"_"+strand
            additionalfields = "0\t" + bedstrand + "\t" + str(refpos+refstart) + "\t" + str(refpos+refstart+oplength+extend) + "\t0,0,0\t" + alignstring
            variantlist.append(bedinterval(chrom=ref, start=refpos+refstart, end=refpos+refstart+oplength+extend, name=variantname, rest=additionalfields ))

        if op == 1:
            refallele = "*"
            queryallele = queryseq[querycurrentoffset:querycurrentoffset+oplength]
            if widen is True:
                extend = 0
                while querycurrentoffset + extend < querylength and refcurrentoffset + extend < reflength and refseq[refcurrentoffset + extend] == queryseq[querycurrentoffset + extend]:
                    queryallele = queryallele + refseq[refcurrentoffset + extend]
                    if refallele == "*":
                        refallele = refseq[refcurrentoffset + extend]
                    else:
                        refallele = refallele + refseq[refcurrentoffset + extend]
                    extend = extend + 1
            if strand == 'F':
                querycoordinate = querystart + querycurrentoffset
            else:
                querycoordinate = queryend - querycurrentoffset
            variantname=query+"_"+str(querycoordinate)+"_"+refallele+"_"+queryallele+"_"+strand
            additionalfields = "0\t" + bedstrand + "\t" + str(refpos+refstart) + "\t" + str(refpos+refstart+extend) + "\t0,0,0\t" + alignstring
            variantlist.append(bedinterval(chrom=ref, start=refpos+refstart, end=refpos+refstart+extend, name=variantname, rest=additionalfields ))


        # advance current positions: cases where reference coord advances (MDN=X):
        if op in [0, 2, 3, 7, 8]:
            refcurrentoffset = refcurrentoffset + oplength
        # cases where query coord advances (MI=X)
        if op in [0, 1, 7, 8]:
            querycurrentoffset = querycurrentoffset + oplength

        alignopindex = alignopindex + 1

    return variantlist

