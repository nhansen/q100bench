import re
import pysam
import pybedtools
import logging
from collections import namedtuple
from pathlib import Path
from q100bench import seqparse
from q100bench import phasing
from q100bench import bedtoolslib
from q100bench import output

# create namedtuple for bed intervals:
varianttuple = namedtuple('varianttuple', ['chrom', 'start', 'end', 'name', 'vartype', 'excluded']) 

logger = logging.getLogger(__name__)

# query start is always less than query end regardless of strand. query left always corresponds to ref start, 
# and so will be greater than query right when strand is reversed--all four are 1-based
def write_bedfiles(bamobj, pafaligns, refobj, queryobj, hetsites, testmatbed, testpatbed, truthbed, hetallelebed, excludedbedobj, args):

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
                query, querystart, queryend, ref, refstart, refend, strand = retrieve_align_data(align)
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
        # mark variants that are in excluded regions:
        variants = exclude_variants(variants, excludedbedobj)
    else:
        for pafdict in pafaligns:
            # all start/endpoints are 1-based
            query = pafdict['query']
            querystart = pafdict['querystart']
            queryend = pafdict['queryend']
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
            querycoveredstring += query + "\t" + str(queryleft - 1) + "\t" + str(queryright) + "\t" + refnamestring + "\n"
            refcoveredstring += ref + "\t" + str(refstart - 1) + "\t" + str(refend) + "\t" + querynamestring + "\n"

    refcoveredbed = pybedtools.BedTool(refcoveredstring, from_string = True)
    querycoveredbed = pybedtools.BedTool(querycoveredstring, from_string = True)

    phap1 = re.compile(r'.*MAT.*')
    phap2 = re.compile(r'.*PAT.*')
    with open(testmatbed, "w") as tmb:
        for testint in sorted(querycoveredbed, key=lambda h: (h.chrom, h.start, h.stop)):
            if phap1.match(testint.name):
                tmb.write(testint.chrom + "\t" + str(testint.start) + "\t" + str(testint.end) + "\t" + testint.name + "\n")

    with open(testpatbed, "w") as tpb:
        for testint in sorted(querycoveredbed, key=lambda h: (h.chrom, h.start, h.stop)):
            if phap2.match(testint.name):
                tpb.write(testint.chrom + "\t" + str(testint.start) + "\t" + str(testint.end) + "\t" + testint.name + "\n")

    with open(truthbed, "w") as rb:
        for truthint in sorted(refcoveredbed, key=lambda h: (h.chrom, h.start, h.stop)):
            rb.write(truthint.chrom + "\t" + str(truthint.start) + "\t" + str(truthint.end) + "\t" + truthint.name + "\n")

    phasing.write_hetallele_bed(hetsitealleles, hetallelebed)

    if user_variantfile is not None:
        with open(user_variantfile, "r") as vh:
            variantline = vh.readline()
            while variantline:
                variantline = variantline.rstrip()
                [chrom, start, end, name] = variantline.split("\t")
                variants.append(varianttuple(chrom=chrom, start=int(start), end=int(end), name=name))
                variantline = vh.readline()

    return [refcoveredbed, querycoveredbed, variants, hetsitealleles]

def retrieve_align_data(align)->list:
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
def align_variants(align, queryobj, query:str, querystart:int, queryend:int, refobj, ref:str, refstart:int, refend:int, strand:str, chromhetsites={}, hetsitealleles={}, widen=True)->list:

    # coordinates are all one-based, with start at beginning of *original* sequence (not left end of the alignment)
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
    refalignlength = refend - refstart + 1
    queryalignlength = queryend - querystart + 1
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
                    #additionalfields = "0\t" + bedstrand + "\t" + str(refpos+refstart-1) + "\t" + str(refpos+refstart) + "\t0,0,0\t" + alignstring
                    variantlist.append(varianttuple(chrom=ref, start=refpos+refstart-1, end=refpos+refstart, name=variantname, vartype='SNV', excluded=False))
                query_positions.append(querypos)

        if op in [2, 3]: # deletions
            refpos = refcurrentoffset-1;
            for deloffset in range(oplength):
                query_positions.append(querycurrentoffset)
            refallele = refseq[refcurrentoffset:refcurrentoffset+oplength] # one-based refstart+refcurrentoffset to refstart+refcurrentoffset+oplength-1
            queryallele = "*" # one-based between querystart+querycurrentoffset-1 and querystart+querycurrentoffset if forward strand, queryend-querycurrentoffset+1 and queryend-querycurrentoffset if rev
            extendleft = 0
            extendright = 0
            if widen is True: # n.b. - this will *lower* the righthand coordinate of reverse strand queries by "extendright" and *lower* the lefthand coordinate of ref entries and forward strand queries by "extendleft"
                while querycurrentoffset + extendright < queryalignlength and refcurrentoffset + extendright < refalignlength and refseq[refcurrentoffset + extendright] == queryseq[querycurrentoffset + extendright]:
                    refallele = refallele + queryseq[querycurrentoffset + extendright]
                    if queryallele == "*":
                        queryallele = queryseq[querycurrentoffset + extendright]
                    else:
                        queryallele = queryallele + queryseq[querycurrentoffset + extendright]
                    extendright = extendright + 1
                while querycurrentoffset - 1 - extendleft >= 0 and refcurrentoffset - 1 + oplength - extendleft >= 0 and refseq[refcurrentoffset - 1 + oplength - extendleft] == queryseq[querycurrentoffset - 1 - extendleft]:
                    refallele = queryseq[querycurrentoffset - extendleft - 1] + refallele
                    if queryallele == "*":
                        queryallele = queryseq[querycurrentoffset - extendleft - 1]
                    else:
                        queryallele = queryseq[querycurrentoffset - extendleft - 1] + queryallele
                    extendleft = extendleft + 1
            # querycoordinate values still need to be adjusted due to widening?
            if strand == 'F':
                querycoordinate = querystart + querycurrentoffset - extendleft
            else:
                querycoordinate = queryend - querycurrentoffset - extendright

            # check neighboring bases for Ns, which can create misleading/possibly wrong variants
            # (looks as if I might be looking only at left and right bases within the query alleles here--need to check)
            [queryleftbase, queryrightbase] = ["", ""]
            if querycurrentoffset > 0:
                queryleftbase = queryseq[querycurrentoffset-1]
            if querycurrentoffset+extendright < queryalignlength:
                queryrightbase = queryseq[querycurrentoffset+extendright]
            querysurroundingseq = queryleftbase + queryrightbase
        
            if not (matchns.match(queryallele) or matchns.match(refallele) or matchns.match(querysurroundingseq)):
                variantname=query+"_"+str(querycoordinate+1)+"_"+refallele+"_"+queryallele+"_"+strand # positions of insertions are positions to the left of first inserted base
                logger.debug("Variant with pos " + ref + ":" + str(refpos+refstart-extendleft) + "-" + str(refpos+refstart+oplength+extendright) + " name " + variantname + " and queryallele " + queryallele + " and refallele " + refallele + " has query surrounding seq " + querysurroundingseq)
                variantlist.append(varianttuple(chrom=ref, start=refpos+refstart-extendleft, end=refpos+refstart+oplength+extendright, name=variantname, vartype='INDEL', excluded=False ))
            else:
                logger.debug("Variant with name " + variantname + " and queryallele " + queryallele + " and refallele " + refallele + " has query surrounding seq " + querysurroundingseq + " was excluded")

        if op == 1: # insertion
            refpos = refcurrentoffset-1;
            refallele = "*"
            queryallele = queryseq[querycurrentoffset:querycurrentoffset+oplength] # one-based querystart+querycurrentoffset to querystart+querycurrentoffset+oplength-1 if forward strand, queryend-querycurrentoffset to queryend-querycurrentoffset-oplength+1 if reverse
            #query_positions.append(querycurrentoffset)
            extendright = 0
            extendleft = 0
            if widen is True: # n.b. - this will *lower* the righthand coordinate of reverse strand queries by "extendright"
                while querycurrentoffset + extendright < queryalignlength and refcurrentoffset + extendright < refalignlength and refseq[refcurrentoffset + extendright] == queryseq[querycurrentoffset + extendright]:
                    queryallele = queryallele + refseq[refcurrentoffset + extendright]
                    if refallele == "*":
                        refallele = refseq[refcurrentoffset + extendright]
                    else:
                        refallele = refallele + refseq[refcurrentoffset + extendright]
                    extendright = extendright + 1
                while querycurrentoffset - 1 + oplength - extendleft >= 0 and refcurrentoffset - 1 - extendleft >= 0 and refseq[refcurrentoffset - 1 - extendleft] == queryseq[querycurrentoffset - 1 + oplength - extendleft]:
                    queryallele = refseq[refcurrentoffset - extendleft - 1] + queryallele
                    if refallele == "*":
                        refallele = refseq[refcurrentoffset - extendleft - 1]
                    else:
                        refallele = refseq[refcurrentoffset - extendleft - 1] + refallele
                    extendleft = extendleft + 1
            if strand == 'F':
                querycoordinate = querystart + querycurrentoffset - extendleft
                querycoordend = querycoordinate + oplength - 1 # this is potentially off by one and could be a bug (see its use below)
            else:
                querycoordinate = queryend - querycurrentoffset - extendleft
                querycoordend = querycoordinate - oplength - 1 # this is potentially off by one and could be a bug (see its use below)

            # check neighboring bases for Ns, which can create misleading/possibly wrong variants
            [refleftbase, refrightbase] = ["", ""]
            if refcurrentoffset-extendleft > 0:
                refleftbase = refseq[refcurrentoffset-1-extendleft]
            else:
                refleftbase = ''
            if refcurrentoffset+extendright+1 < refalignlength:
                refrightbase = refseq[refcurrentoffset+extendright+1]
            else:
                refrightbase = ''
            refsurroundingseq = refleftbase + refrightbase
            
            if not (matchns.match(queryallele) or matchns.match(refallele) or matchns.match(refsurroundingseq)):
                variantname=query+"_"+str(querycoordinate)+"_"+refallele+"_"+queryallele+"_"+strand
                #additionalfields = "0\t" + bedstrand + "\t" + str(refpos+refstart) + "\t" + str(refpos+refstart+extendright) + "\t0,0,0\t" + alignstring
                logger.debug("Variant with pos " + ref + ":" + str(refpos+refstart-extendleft) + "-" + str(refpos+refstart+extendright) + " name " + variantname + " and queryallele " + queryallele + " and refallele " + refallele + " has ref surrounding seq " + refsurroundingseq)
                variantlist.append(varianttuple(chrom=ref, start=refpos+refstart-extendleft, end=refpos+refstart+extendright, name=variantname, vartype='INDEL', excluded=False ))
            else:
                logger.debug("Variant with name " + variantname + " and queryallele " + queryallele + " and refallele " + refallele + " has ref surrounding seq " + refsurroundingseq + " was excluded")

        # advance current positions: cases where reference coord advances (MDN=X):
        if op in [0, 2, 3, 7, 8]:
            refcurrentoffset = refcurrentoffset + oplength
        # cases where query coord advances (MI=X)
        if op in [0, 1, 7, 8]:
            querycurrentoffset = querycurrentoffset + oplength

        alignopindex = alignopindex + 1

    numquerypositions = len(query_positions)
    refseqlength = len(refseq)

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

    return variantlist

def exclude_variants(variants:list, excludedregionsobj:pybedtools.BedTool)->list:
    # create bedintervals for variants:
    variantbedstring = ""
    for variant in variants:
        chrom = variant.chrom
        start = variant.start
        end = variant.end
        name = variant.name

        variantbedstring = variantbedstring + chrom + "\t" + str(start) + "\t" + str(end) + "\t" + name + "\n"

    variantbedobj = pybedtools.BedTool(variantbedstring, from_string=True)
    excludedvariants = bedtoolslib.intersectintervals(variantbedobj, excludedregionsobj, wa=True)
    excludeddict = {}
    for excludedvariant in excludedvariants:
        excludeddict[excludedvariant.name] = True

    newvariants = []
    for variant in variants:
        if variant.name in excludeddict.keys():
            newvariants.append(varianttuple(chrom=variant.chrom, start=variant.start, end=variant.end, name=variant.name, vartype=variant.vartype, excluded = True))
        else:
            newvariants.append(variant)

    return newvariants

def read_paf_aligns(paffile:str, mintargetlength=0)->list:

    alignlist = []
    
    paf = Path(paffile)
    if not paf.is_file():
        logger.critical("PAF file " + paffile + " must exist and be readable")
        exit(1)

    with open(paffile, "r") as pfh:
        alignline = pfh.readline()
        while alignline:
            alignline = alignline.rstrip()
            fields = alignline.split("\t")
            if len(fields) >= 12:
                [query, querylength, querystartzb, queryend, strand, target, targetlength, targetstartzb, targetend, resmatches, blocklength, mapqual] = fields[0:12]
            else:
                logger.critical("Input paf-file has fewer than 12 tab-delimited columns. Unable to process.")
                exit(1)

            #if len(fields) > 12 and len(list(filter(lambda x:'tp:A:S' in x, fields[12:]))) > 0:
                #continue

            querystart = int(querystartzb) + 1
            queryend = int(queryend)
            queryalignlength = queryend - int(querystartzb)
            targetstart = int(targetstartzb) + 1
            targetend = int(targetend)
            targetalignlength = targetend - int(targetstartzb)

            if targetalignlength >= mintargetlength:
                if strand == "+":
                    alignlist.append({'query':query, 'querylength':int(querylength), 'querystart':querystart, 'queryend':queryend, 'queryalignlength':queryalignlength, 'strand':strand, 'target':target, 'targetlength':int(targetlength), 'targetstart':targetstart, 'targetend':targetend, 'targetalignlength':targetalignlength})
                else:
                    alignlist.append({'query':query, 'querylength':int(querylength), 'queryend':querystart, 'querystart':queryend, 'queryalignlength':queryalignlength, 'strand':strand, 'target':target, 'targetlength':int(targetlength), 'targetstart':targetstart, 'targetend':targetend, 'targetalignlength':targetalignlength})
            alignline = pfh.readline()

    return alignlist

def read_bam_aligns(bamobj, mintargetlength=0)->list:

    alignlist = []

    refentries = bamobj.references
    reflengths = bamobj.lengths
    reflengthdict = {}
    for tid in range(len(refentries)):
        reflengthdict[refentries[tid]] = reflengths[tid]

    for align in bamobj.fetch():
        if align.is_secondary:
            continue
        if align.reference_length >= mintargetlength:
            query, querystart, queryend, ref, refstart, refend, strand = retrieve_align_data(align)
            querylength = align.query_length
            targetlength = reflengthdict[ref]
            if strand == "F":
                queryleft = querystart
                queryright = queryend
                alignlist.append({'query':query, 'querylength':querylength, 'querystart':querystart, 'queryend':queryend, 'queryalignlength':queryend-querystart+1, 'strand':'+', 'target':ref, 'targetlength':targetlength, 'targetstart':refstart, 'targetend':refend, 'targetalignlength':refend-refstart+1})
            else:
                queryleft = queryend
                queryright = querystart
                alignlist.append({'query':query, 'querylength':querylength, 'queryend':querystart, 'querystart':queryend, 'queryalignlength':queryend-querystart+1, 'strand':'-', 'target':ref, 'targetlength':targetlength, 'targetstart':refstart, 'targetend':refend, 'targetalignlength':refend-refstart+1})

    return alignlist

# this routine assumes query start > query end for reverse strand alignments
def assess_overall_structure(aligndata:list, refobj, queryobj, outputfiles, bedobjects, benchmark_stats, args):

    # create a "clustercoverage" dictionary that will contain each chromosome's number of clusters and bases covered
    benchmark_stats["clustercoverage"] = {}
    benchmark_stats["alignclusters"] = {}
    benchmark_stats["numnonexcludedbases"] = {}

    # maximum separating distance along target to include in one cluster of alignments
    maxdistance = args.maxclusterdistance

    # create a directory for alignment line plots for each chromosome:
    output.create_output_directory(outputfiles["alignplotdir"])
    alignplotprefix = outputfiles["alignplotprefix"]

    # calculate number of non-excluded bases for each ref entry:
    for ref in refobj.references:
        benchmark_stats["numnonexcludedbases"][ref] = refobj.get_reference_length(ref)

    # since allexcludedregions bed object is merged,can subtract out excluded bases interval by interval
    for interval in bedobjects["allexcludedregions"]:
        ref = interval.chrom
        benchmark_stats["numnonexcludedbases"][ref] = benchmark_stats["numnonexcludedbases"][ref] - len(interval)

    aligndict = {}
    for align in aligndata:
        refentry = align["target"]
        if refentry not in aligndict:
            aligndict[refentry] = [align]
        else:
            aligndict[refentry].append(align)

    # assess each benchmark entry, one by one
    for refentry in sorted(aligndict.keys()):
        refnelength = benchmark_stats["numnonexcludedbases"][refentry]
        refalignclusters = []
        # sort alignments from longest (along the benchmark) to shortest:
        aligndict[refentry].sort(reverse=True, key=lambda align: align["targetalignlength"])
        # calculate slope as query diff over ref diff, and cluster alignments within the same query/target band:
        for refalign in aligndict[refentry]:
            alignrefstart = refalign['targetstart']
            alignrefend = refalign['targetend']
            alignquery = refalign['query']
            alignquerystart = refalign['querystart']
            alignqueryend = refalign['queryend']
            add_align_to_clusters(refalign, refalignclusters, maxdistance)

        # split clusters that are separated along the target by more than maxdistance:
        refalignclusters = split_disjoint_clusters(refalignclusters, maxdistance)

        for cluster in sorted(refalignclusters, key=lambda c: c["aligns"][0]["targetstart"]):
            clusterquery = cluster["query"] 
            clusterslope = cluster["slope"]
            clusterintercept = cluster["intercept"]
            clusterbedstring = ""
            for align in sorted(cluster["aligns"], key=lambda a: (a["targetstart"], a["targetend"])):
                clusterbedstring = clusterbedstring + refentry + "\t" + str(align["targetstart"]) + "\t" + str(align["targetend"]) + "\t" + clusterquery + "_" + str(align["querystart"]) + "_" + str(align["queryend"]) + "\n"

            bedtool = pybedtools.BedTool(clusterbedstring, from_string = True)
            nonexcludedbedtool = bedtoolslib.intersectintervals(bedtool, bedobjects["allexcludedregions"], v=True)
            cluster["nonexcludedcoveredbases"] = bedtoolslib.bedsum(nonexcludedbedtool)
            if cluster["nonexcludedcoveredbases"] > 0:
                mergednonexcludedbedtool = bedtoolslib.mergeintervals(nonexcludedbedtool)
                cluster["nrnonexcludedcoveredbases"] = bedtoolslib.bedsum(mergednonexcludedbedtool)
            else:
                cluster["nrnonexcludedcoveredbases"] = 0

        # calculate how many clusters needed to cover 95% of ref:
        totalnonexcludedcovered = 0
        clustercount = 0
        lca95 = None
        nca95 = None
        clusterno = 1
        refentrybedstring = ""
        for cluster in sorted(refalignclusters, key=lambda c:c["nrnonexcludedcoveredbases"], reverse = True):
            clusterquery = cluster["query"] 
            for align in sorted(cluster["aligns"], key=lambda a: (a["targetstart"], a["targetend"])):
                if lca95 is None:
                    clustername = "Cluster" + str(clusterno)
                else:
                    clustername = "SmallCluster" + str(clusterno)
                refentrybedstring = refentrybedstring + refentry + "\t" + str(align["targetstart"]) + "\t" + str(align["targetend"]) + "\t" + clusterquery + "_" + str(align["querystart"]) + "_" + str(align["queryend"]) + "_" + clustername + "\n"
            clusterno = clusterno + 1

            if lca95 is None:
                totalnonexcludedcovered = totalnonexcludedcovered + cluster["nrnonexcludedcoveredbases"]
            clustercount = clustercount + 1
            if lca95 is None and totalnonexcludedcovered > 0.95*refnelength:
                lca95 = clustercount
                nca95 = cluster["nrnonexcludedcoveredbases"]

        refentrybedtool = pybedtools.BedTool(refentrybedstring, from_string = True)
        refentrybedtool.saveas(alignplotprefix + "." + refentry + ".clusters.bed")

        # Note: lca95 will be "None" for reference entries not able to be covered
        benchmark_stats["clustercoverage"][refentry] = {"lca95":lca95, "nonexcludedcovered":totalnonexcludedcovered}
        refalignclusters.sort(key=lambda c: c["aligns"][0]["targetstart"])
        benchmark_stats["alignclusters"][refentry] = refalignclusters
    
    return benchmark_stats

def add_align_to_clusters(align:dict, alignclusters:list, maxdistance:int):

    alignstart = align['targetstart']
    alignend = align['targetend']
    alignquery = align['query']
    alignquerystart = align['querystart']
    alignqueryend = align['queryend']
    alignslope = (alignend - alignstart)/(alignqueryend-alignquerystart)
    alignintercept = alignstart - int(alignslope * alignquerystart)

    # try to assign this align to a pre-existing cluster of aligns:
    assigned = False
    for cluster in alignclusters:
        if cluster["query"] != alignquery:
            continue
        clusterquery = cluster["query"]
        clusterslope = cluster["slope"]
        clusterintercept = cluster["intercept"]
        predstart = clusterintercept + clusterslope * alignquerystart
        if abs(predstart - alignstart) <= maxdistance:
            cluster["aligns"].append(align)
            assigned = True
            break

    # create a new cluster if none were appropriate
    if not assigned:
        alignclusters.append({'query':alignquery, 'slope':alignslope, 'intercept':alignintercept, 'aligns':[align]})

    return 0

def split_disjoint_clusters(refalignclusters:list, maxdistance:int)->list:

    spinoffclusters = []
    for cluster in refalignclusters:
        if len(cluster["aligns"])==1:
            continue

        runningalignlist = []
        maxposition = 0
        for align in sorted(cluster["aligns"], key=lambda a: (a["targetstart"], a["targetend"])):
            # do we have too big a gap?
            if len(runningalignlist) > 0 and align["targetstart"] - maxposition > maxdistance:
                spinoffclusters.append({'query':cluster["query"], 'slope':cluster["slope"], 'intercept':cluster["intercept"], 'aligns':runningalignlist})
                runningalignlist = []
                maxposition = 0
            runningalignlist.append(align)
            maxposition = max(maxposition, align["targetend"])
        cluster["aligns"] = runningalignlist

    return refalignclusters + spinoffclusters

