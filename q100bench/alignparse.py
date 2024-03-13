import re
import pysam
import pybedtools
from collections import namedtuple
from q100bench import seqparse
from q100bench import phasing

# create namedtuple for bed intervals:
bedinterval = namedtuple('bedinterval', ['chrom', 'start', 'end', 'name', 'rest']) 

# query start is always less than query end regardless of strand. query left always corresponds to ref start, 
# and so will be greater than query right when strand is reversed--all four are 1-based
def write_bedfiles(bamobj, pafaligns, refobj, queryobj, hetsites, testmatbed, testpatbed, truthbed, variantbed, hetallelebed, args):

    refcoveredstring = ""
    querycoveredstring = ""
    user_variantfile = args.variantfile
    variants = []
    hetsitealleles = {}

    if bamobj is not None:
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
                querycoveredstring += query + "\t" + str(querystart - 1) + "\t" + str(queryend) + "\t" + refnamestring + "\n"
                refcoveredstring += ref + "\t" + str(refstart - 1) + "\t" + str(refend) + "\t" + querynamestring + "\n"
                if user_variantfile is None:
                    variants.extend(align_variants(align, queryobj, query, querystart, queryend, refobj, ref, refstart, refend, strand, hetsites, hetsitealleles, True))
    else:
        for pafdict in pafaligns:
            # all start/endpoints are 1-based
            query = pafdict['query']
            querystart = pafdict['querystart']
            queryend = pafdict['queryend']
            #print(query + ":" + str(querystart) + "-" + str(queryend))
            if querystart < queryend:
                queryleft = querystart
                queryright = queryend
            else:
                queryleft = queryend
                queryright = querystart
            ref = pafdict['target']
            refstart = pafdict['targetstart']
            refend = pafdict['targetend']
            strand = pafdict['strand']

            querynamestring = query + "." + str(querystart) + "." + str(queryend)
            refnamestring = ref + "." + str(refstart) + "." + str(refend) + "." + strand
            #print(querynamestring + ", " + refnamestring)
            querycoveredstring += query + "\t" + str(queryleft - 1) + "\t" + str(queryright) + "\t" + refnamestring + "\n"
            refcoveredstring += ref + "\t" + str(refstart - 1) + "\t" + str(refend) + "\t" + querynamestring + "\n"

    refcoveredbed = pybedtools.BedTool(refcoveredstring, from_string = True)
    querycoveredbed = pybedtools.BedTool(querycoveredstring, from_string = True)

    phap1 = re.compile(r'.*MAT.*')
    phap2 = re.compile(r'.*PAT.*')
    with open(testmatbed, "w") as tmb:
        for testtuple in sorted(querycoveredbed, key=lambda h: (h.chrom, h.start, h.stop)):
            if phap1.match(testtuple.name):
                tmb.write(testtuple.chrom + "\t" + str(testtuple.start) + "\t" + str(testtuple.end) + "\t" + testtuple.name + "\n")

    with open(testpatbed, "w") as tpb:
        for testtuple in sorted(querycoveredbed, key=lambda h: (h.chrom, h.start, h.stop)):
            if phap2.match(testtuple.name):
                tpb.write(testtuple.chrom + "\t" + str(testtuple.start) + "\t" + str(testtuple.end) + "\t" + testtuple.name + "\n")

    with open(truthbed, "w") as rb:
        for truthtuple in sorted(refcoveredbed, key=lambda h: (h.chrom, h.start, h.stop)):
            rb.write(truthtuple.chrom + "\t" + str(truthtuple.start) + "\t" + str(truthtuple.end) + "\t" + truthtuple.name + "\n")

    phasing.write_hetallele_bed(hetsitealleles, hetallelebed)

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
                variants.append(bedinterval(chrom=chrom, start=int(start), end=int(end), name=name, rest=''))
                variantline = vh.readline()

    return [refcoveredbed, querycoveredbed, variants]

def retrieve_align_data(align, args)->list:
    if align.is_reverse:
        strand = 'R'
    else:
        strand = 'F'
    query = align.query_name

    # number of hard clipped bases:
    cigartuples = align.cigartuples
    # will use actual one-based positions, so I don't go crazy, then report BED format zero-based half open later
    # pysam's align.query_alignment_start is the 0-based coordinate within the uncomplemented hard clipped query sequence, so here we add hardclipping from the 
    # appropriate end to get coordinates within the entire sequence (so low coordinates are towards the start of the original query sequence, not the left end
    # of the alignment)
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
def align_variants(align, queryobj, query:str, querystart:int, queryend:int, refobj, ref:str, refstart:int, refend:int, strand:str, chromhetsites:dict, hetsitealleles:dict, widen=True)->list:

    # coordinates are all one-based, with start at beginning of *original* sequence (not left end of the alignment)
    #print("In align_variants with " + query + ":" + str(querystart) + "-" + str(queryend) + " " + ref + ":" + str(refstart) + "-" + str(refend) + "/" + strand)
    variantlist = []
    coveredregionlist = []
    homozygousregionlist = []

    # make an array of query positions for each ref position:
    query_positions = []

    if queryobj is None:
        queryseq = align.query_alignment_sequence
    else:
        queryseq = queryobj.fetch(reference=query, start=querystart-1, end=queryend).upper()

    refseq = refobj.fetch(reference=ref, start=refstart-1, end=refend).upper()

    alignstring = ref + ":" + str(refstart-1) + "-" + str(refend) + ";" + query + ":" + str(querystart-1) + "-" + str(queryend)
    strandsign = 1
    bedstrand = '+'
    if strand == 'R':
        if queryobj is not None:
            queryseq = seqparse.revcomp(queryseq)
        strandsign = -1
        bedstrand = '-'

    alignops = align.cigartuples

    # first position will begin at left-most ref/query base of the alignment (regardless of strand)
    refcurrentoffset = 0
    querycurrentoffset = 0

    alignopindex = 0
    reflength = refend - refstart + 1
    querylength = queryend - querystart + 1
    matchns = re.compile(".*[nN].*")

    while refcurrentoffset <= refend-refstart and alignopindex < len(alignops): # traverse the alignment operator by operator
        alignop = alignops[alignopindex]
        op = alignop[0]
        oplength = alignop[1]

        if op in [0, 7, 8]: # MX= find SNV and MNVs
            for blockoffset in range(oplength):
                refpos = refcurrentoffset + blockoffset # this is distance from left-most base of the alignment
                querypos = querycurrentoffset + blockoffset # this is distance from left-most base of the alignment
                if strand == 'F':
                    querycoordinate = querypos + querystart
                else:
                    querycoordinate = queryend - querypos
                if refseq[refpos] != queryseq[querypos] and refseq[refpos] != "N" and queryseq[querypos] != "N":
                    variantname=query+"_"+str(querycoordinate)+"_"+refseq[refpos]+"_"+queryseq[querypos]+"_"+strand # query's 1-based position, ref base, query base (comp if rev strand), strand
                    additionalfields = "0\t" + bedstrand + "\t" + str(refpos+refstart-1) + "\t" + str(refpos+refstart) + "\t0,0,0\t" + alignstring
                    variantlist.append(bedinterval(chrom=ref, start=refpos+refstart-1, end=refpos+refstart, name=variantname, rest=additionalfields))
                query_positions.append(querypos)

        if op in [2, 3]: # deletions
            refpos = refcurrentoffset-1;
            for deloffset in range(oplength):
                query_positions.append(querycurrentoffset)
            refallele = refseq[refcurrentoffset:refcurrentoffset+oplength] # one-based refstart+refcurrentoffset to refstart+refcurrentoffset+oplength-1
            queryallele = "*" # one-based between querystart+querycurrentoffset-1 and querystart+querycurrentoffset if forward strand, queryend-querycurrentoffset+1 and queryend-querycurrentoffset if rev
            if widen is True: # n.b. - this will *lower* the righthand coordinate of reverse strand queries by "extend"
                extend = 0
                while querycurrentoffset + extend < querylength and refcurrentoffset + extend < reflength and refseq[refcurrentoffset + extend] == queryseq[querycurrentoffset + extend]:
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

            # check neighboring bases for Ns, which can create misleading/possibly wrong variants
            [queryleftbase, queryrightbase] = ["", ""]
            if querycurrentoffset > 0:
                queryleftbase = queryseq[querycurrentoffset-1]
            if querycurrentoffset+extend < querylength:
                queryrightbase = queryseq[querycurrentoffset+extend]
            querysurroundingseq = queryleftbase + queryrightbase
        
            if not (matchns.match(queryallele) or matchns.match(refallele) or matchns.match(querysurroundingseq)):
                variantname=query+"_"+str(querycoordinate+1)+"_"+refallele+"_"+queryallele+"_"+strand # positions of insertions are positions to the left of first inserted base
                additionalfields = "0\t" + bedstrand + "\t" + str(refpos+refstart) + "\t" + str(refpos+refstart+oplength+extend) + "\t0,0,0\t" + alignstring
                variantlist.append(bedinterval(chrom=ref, start=refpos+refstart, end=refpos+refstart+oplength+extend, name=variantname, rest=additionalfields ))

        if op == 1: # insertion
            refpos = refcurrentoffset-1;
            refallele = "*"
            queryallele = queryseq[querycurrentoffset:querycurrentoffset+oplength] # one-based querystart+querycurrentoffset to querystart+querycurrentoffset+oplength-1 if forward strand, queryend-querycurrentoffset to queryend-querycurrentoffset-oplength+1 if reverse
            #query_positions.append(querycurrentoffset)
            if widen is True: # n.b. - this will *lower* the righthand coordinate of reverse strand queries by "extend"
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
                querycoordend = querycoordinate + oplength - 1 # this is potentially off by one and could be a bug (see its use below)
            else:
                querycoordinate = queryend - querycurrentoffset
                querycoordend = querycoordinate - oplength - 1 # this is potentially off by one and could be a bug (see its use below)

            # check neighboring bases for Ns, which can create misleading/possibly wrong variants
            [refleftbase, refrightbase] = ["", ""]
            if refcurrentoffset > 0:
                refleftbase = refseq[refcurrentoffset-1]
            if refcurrentoffset+extend+1 < reflength:
                refrightbase = refseq[refcurrentoffset+extend+1]
            refsurroundingseq = refleftbase + refrightbase
            
            if not (matchns.match(queryallele) or matchns.match(refallele) or matchns.match(refsurroundingseq)):
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

    numquerypositions = len(query_positions)
    refseqlength = len(refseq)
    #print(str(numquerypositions) + " positions for " + str(refseqlength) + " bases")

    # use query positions to assess het alleles/covered regions:
    if ref in chromhetsites:
        desiredhets = chromhetsites[ref]
        for het in desiredhets:
            hetname = het.name
            namefields = hetname.split("_")
            hetpos = int(namefields[-4])
            refallele = namefields[-3]
            altallele = namefields[-2]

            hetstart = het.start
            hetend = het.end
            if hetstart < refstart:
                continue
            if hetend >= refend:
                break
            querystartoffset = query_positions[hetstart-refstart]
            if hetend - refstart + 2 >= len(query_positions):
                queryendoffset = len(queryseq)
            else:
                queryendoffset = query_positions[hetend-refstart+2]
            queryallele = '*'
            if querystartoffset + 1 < queryendoffset - 1:
                queryallele = queryseq[querystartoffset+1:queryendoffset-1]
            alleletype = 'NEITHER'
            if queryallele == refallele:
                alleletype = 'SAME'
            elif queryallele == altallele:
                alleletype = 'ALT'

            querystartcoord = querystart + querystartoffset
            queryendcoord = querystart + queryendoffset - 2

            hetsitealleles[hetname] = {'name':hetname, 'ref':ref, 'refstart':hetstart, 'refend':hetend, 'allele':queryallele, 'query':query, 'start':querystartcoord, 'end':queryendcoord, 'chrom':query}
            #print(hetname + "\t" + ref + "\t" + str(hetstart) + "\t" + str(hetend) + "\t" + query + "\t" + str(querystartcoord) + "\t" + str(queryendcoord) + "\t" + queryallele + "\t" + alleletype)

    return variantlist

def read_paf_aligns(paffile:str, mintargetlength)->list:

    alignlist = []
    with open(paffile, "r") as pfh:
        alignline = pfh.readline()
        while alignline:
            alignline = alignline.rstrip()
            fields = alignline.split("\t")
            if len(fields) >= 12:
                [query, querylength, querystartzb, queryend, strand, target, targetlength, targetstartzb, targetend, resmatches, blocklength, mapqual] = fields[0:12]
            else:
                sys.print("Input paf-file has fewer than 12 tab-delimited columns. Unable to process.")
                exit(1)
            if int(targetend) - int(targetstartzb) >= mintargetlength:
                print("Keeping paf entry " + target + ":" + targetstartzb + "-" + targetend)
                if strand == "+":
                    alignlist.append({'query':query, 'querylength':int(querylength), 'querystart':int(querystartzb)+1, 'queryend':int(queryend), 'strand':strand, 'target':target, 'targetlength':int(targetlength), 'targetstart':int(targetstartzb)+1, 'targetend':int(targetend)})
                else:
                    alignlist.append({'query':query, 'querylength':int(querylength), 'queryend':int(querystartzb)+1, 'querystart':int(queryend), 'strand':strand, 'target':target, 'targetlength':int(targetlength), 'targetstart':int(targetstartzb)+1, 'targetend':int(targetend)})
            alignline = pfh.readline()

    return alignlist
