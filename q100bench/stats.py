import re
import math

def write_general_assembly_stats(generalstatspath, refobj, queryobj, contiglist:list, gaplist:list, args)->dict:

    bmstats = {}

    # total bases on benchmark haplotypes (for NG/LG calculation):
    hap1totalbases = 0
    hap2totalbases = 0
    phap1 = re.compile(r'.*MAT.*')
    phap2 = re.compile(r'.*PAT.*')
    for refcontig in refobj.references:
        hap1match = phap1.match(refcontig)
        hap2match = phap2.match(refcontig)
        if hap1match:
            hap1totalbases = hap1totalbases + refobj.get_reference_length(refcontig)
        if hap2match:
            hap2totalbases = hap2totalbases + refobj.get_reference_length(refcontig)

    bmstats['hap1totalbases'] = hap1totalbases # this is the total number of bases in the MATERNAL benchmark chroms
    bmstats['hap2totalbases'] = hap2totalbases # this is the total number of bases in the PATERNAL benchmark chroms
    bmstats['totalbases'] = hap1totalbases + hap2totalbases

    mincontiglength = args.mincontiglength
    with open(generalstatspath, "w") as gsfh:
        numscaffolds = queryobj.nreferences
        totalscaffoldbases = sum(queryobj.lengths)

        cumscaffbases = 0
        numscaffs = 0
        scaffold_n50 = None
        scaffold_hap1_ng50 = None
        scaffold_hap2_ng50 = None
        scaffold_l50 = None
        scaffold_hap1_lg50 = None
        scaffold_hap2_lg50 = None
        scaffold_lengths = queryobj.lengths

        scaffold_lengths.sort(reverse=True)
        for scafflength in scaffold_lengths:
            cumscaffbases = cumscaffbases + scafflength
            numscaffs = numscaffs + 1
            if cumscaffbases > 0.5*totalscaffoldbases and scaffold_n50 is None:
                scaffold_n50 = scafflength
                scaffold_l50 = numscaffs
            if cumscaffbases > 0.5*hap1totalbases and scaffold_hap1_ng50 is None:
                scaffold_hap1_ng50 = scafflength
                scaffold_hap1_lg50 = numscaffs
            if cumscaffbases > 0.5*hap2totalbases and scaffold_hap2_ng50 is None:
                scaffold_hap2_ng50 = scafflength
                scaffold_hap2_lg50 = numscaffs

        totalns = 0
        for gap in gaplist:
            gaplength = gap.end - gap.start
            totalns = totalns + gaplength

        numcontigs = len(contiglist)
        contig_n50 = None
        contig_hap1_ng50 = None
        contig_hap2_ng50 = None
        contig_l50 = None
        contig_hap1_lg50 = None
        contig_hap1_lg50 = None
        numlargecontigs = 0
        sizelist = []
        totalsize = 0
        for contig in contiglist:
            contiglength = contig.end - contig.start
            if contiglength >= mincontiglength:
                numlargecontigs = numlargecontigs + 1
                sizelist.append(contiglength)
                totalsize = totalsize + contiglength

        sizelist.sort(reverse=True)
        cumcontigbases = 0
        numcontigs = 0
        for contigsize in sizelist:
            cumcontigbases = cumcontigbases + contigsize
            numcontigs = numcontigs + 1
            if cumcontigbases > 0.5*totalsize and contig_n50 is None:
                contig_n50 = contigsize
                contig_l50 = numcontigs
            if cumcontigbases > 0.5*hap1totalbases and contig_hap1_ng50 is None:
                contig_hap1_ng50 = contigsize
                contig_hap1_lg50 = numcontigs
            if cumcontigbases > 0.5*hap2totalbases and contig_hap2_ng50 is None:
                contig_hap2_ng50 = contigsize
                contig_hap2_lg50 = numcontigs
        
        bmstats['totalns'] = totalns

        bmstats['numscaffolds'] = numscaffolds
        bmstats['totalscaffoldbases'] = totalscaffoldbases
        bmstats['scaffoldn50'] = scaffold_n50
        bmstats['scaffoldhap1ng50'] = scaffold_hap1_ng50
        bmstats['scaffoldhap2ng50'] = scaffold_hap2_ng50
        bmstats['scaffoldl50'] = scaffold_l50
        bmstats['scaffoldhap1lg50'] = scaffold_hap1_lg50
        bmstats['scaffoldhap2lg50'] = scaffold_hap2_lg50

        bmstats['numcontigs'] = numcontigs
        bmstats['numlargecontigs'] = numcontigs
        bmstats['totallargecontigbases'] = totalsize
        bmstats['contign50'] = contig_n50
        bmstats['contighap1ng50'] = contig_hap1_ng50
        bmstats['contighap2ng50'] = contig_hap2_ng50
        bmstats['contigl50'] = contig_l50
        bmstats['contighap1lg50'] = contig_hap1_lg50
        bmstats['contighap2lg50'] = contig_hap2_lg50

        gsfh.write("Scaffolds:\n\n")
        gsfh.write("Number of scaffolds: " + str(numscaffolds) + "\n")
        gsfh.write("Total bases in scaffolds: " + str(totalscaffoldbases) + "\n")
        gsfh.write("Ns per kb in scaffolds: " + str(int(totalns*1000/totalscaffoldbases)) + "\n")
        gsfh.write("Scaffold N50: " + str(round(scaffold_n50/1000000, 3)) + "Mb\n")
        gsfh.write("Scaffold NG50 (MAT): " + str(round(scaffold_hap1_ng50/1000000, 3)) + "Mb\n")
        gsfh.write("Scaffold L50: " + str(scaffold_l50) + "\n")
        gsfh.write("Scaffold LG50 (MAT): " + str(scaffold_hap1_lg50) + "\n")
        gsfh.write("\nContigs:\n\n")
        gsfh.write("Number of contigs: " + str(numcontigs) + "\n")
        gsfh.write("Number of contigs >= " + str(mincontiglength) + ": " + str(numlargecontigs) + "\n")
        gsfh.write("Total bases in contigs >= " + str(mincontiglength) + ": " + str(totalsize) + "\n")
        gsfh.write("Contig N50: " + str(round(contig_n50/1000000, 3)) + "Mb\n")
        gsfh.write("Contig NG50 (MAT): " + str(round(contig_hap1_ng50/1000000, 3)) + "Mb\n")
        gsfh.write("Contig L50: " + str(contig_l50) + "\n")
        gsfh.write("Contig LG50 (MAT): " + str(contig_hap1_lg50) + "\n")
        gsfh.write("\n")

    return bmstats

def write_aligned_stats(refobj, queryobj, truthints, mergedtruthints, mergedtestmatints, mergedtestpatints, bedfiles:dict, bmstats:dict, args)->dict:

    generalstatsfile = bedfiles["generalstatsfile"]

    phap1 = re.compile(r'.*MAT.*')
    phap2 = re.compile(r'.*PAT.*')
    totalbenchcovered = 0

    totalrefaligned = 0
    numberrefaligns = 0

    hap1totalbases = bmstats['hap1totalbases'] # total MATERNAL bases in benchmark

    hap1_nga50 = None
    hap1_lga50 = None
    for truthint in sorted(mergedtruthints, key=lambda h: int(h[2]) - int(h[1]), reverse=True):
        alignlength = int(truthint[2]) - int(truthint[1])
        totalrefaligned = totalrefaligned + alignlength
        numberrefaligns = numberrefaligns + 1

        if totalrefaligned >= 0.5*hap1totalbases and hap1_nga50 is None:
            hap1_nga50 = alignlength
            hap1_lga50 = numberrefaligns

    matbenchcovered = 0
    patbenchcovered = 0
    for truthint in mergedtruthints:
        [chrom, start, end, name] = truthint
        start = int(start)
        end = int(end)
        totalbenchcovered = totalbenchcovered + end - start
        if phap1.match(chrom):
            matbenchcovered = matbenchcovered + end - start
        if phap2.match(chrom):
            patbenchcovered = patbenchcovered + end - start

    totaltestmatcovered = 0
    totaltestpatcovered = 0
    for testint in mergedtestmatints:
        [chrom, start, end, name] = testint
        totaltestmatcovered = totaltestmatcovered + int(end) - int(start)

    for testint in mergedtestpatints:
        [chrom, start, end, name] = testint
        totaltestpatcovered = totaltestpatcovered + int(end) - int(start)

    bmstats["testmattotalcovered"] = totaltestmatcovered
    bmstats["testpattotalcovered"] = totaltestpatcovered
    bmstats["benchtotalcovered"] = totalbenchcovered

    with open(generalstatsfile, "a") as gsfh:
        gsfh.write("\nAligned contig bases:\n\n")
        gsfh.write("NGA50 (for MATERNAL benchmark haplotype): " + str(round(hap1_nga50/1000000, 3)) + "Mb\n")
        gsfh.write("LGA50 (for MATERNAL benchmark haplotype): " + str(hap1_lga50) + "\n")
        perctestmatcovered = int(totaltestmatcovered * 1000 / bmstats['totallargecontigbases'] + 0.5) / 10
        perctestpatcovered = int(totaltestpatcovered * 1000 / bmstats['totallargecontigbases'] + 0.5) / 10
        if bmstats['hap1totalbases'] > 0 or bmstats['hap2totalbases'] > 0:
            percbenchcovered = int(1000 * totalbenchcovered / (bmstats['hap1totalbases'] + bmstats['hap2totalbases']) + 0.5) / 10
        else:
            percbenchcovered = 'NA'
        if bmstats['hap1totalbases'] > 0:
            percmatbenchcovered = int(1000 * matbenchcovered / bmstats['hap1totalbases'] + 0.5) / 10
        else:
            percmatbenchcovered = 'NA'
        if bmstats['hap2totalbases'] > 0:
            percpatbenchcovered = int(1000 * patbenchcovered / bmstats['hap2totalbases'] + 0.5) / 10
        else:
            percpatbenchcovered = 'NA'
        gsfh.write("Total " + args.assembly + " covered with aligns to MAT chromosomes: " + str(totaltestmatcovered) + "/" + str(bmstats['totallargecontigbases']) + " (" + str(perctestmatcovered) + "% of total bases in contigs >= " + str(args.mincontiglength) + " bases)" + "\n")
        gsfh.write("Total " + args.assembly + " covered with aligns to PAT chromosomes: " + str(totaltestpatcovered) + "/" + str(bmstats['totallargecontigbases']) + " (" + str(perctestpatcovered) + "% of total bases in contigs >= " + str(args.mincontiglength) + " bases)" + "\n")
        gsfh.write("Total " + args.benchmark + " covered: " + str(totalbenchcovered) + "/" + str(bmstats['totalbases']) + " (" + str(percbenchcovered) + "% of all MAT and PAT bases in benchmark)" + "\n" )
        gsfh.write("MAT " + args.benchmark + " covered: " + str(matbenchcovered) + "/" + str(bmstats['hap1totalbases']) + " (" + str(percmatbenchcovered) + "% of MAT bases in benchmark)" + "\n")
        gsfh.write("PAT " + args.benchmark + " covered: " + str(patbenchcovered) + "/" + str(bmstats['hap2totalbases']) + " (" + str(percpatbenchcovered) + "% of PAT bases in benchmark)" + "\n")

    return bmstats

def write_qv_stats(benchmark_stats:dict, variantfile:str, bedfiles:dict, args):

    totalerrors = 0
    totalphasingerrors = 0
    totalconsensuserrors = 0

    with open(variantfile, "r") as vfh:
        errorline = vfh.readline()
        while errorline:
            errorline = errorline.rstrip()
            [chrom, start, end, name, errortype] = errorline.split("\t")
            namefields = name.split("_")
            refallele = namefields[-2]
            altallele = namefields[-1]

            totalerrors = totalerrors + 1
            if errortype == "PHASING":
                totalphasingerrors = totalphasingerrors + 1
            elif errortype == "CONSENSUS":
                totalconsensuserrors = totalconsensuserrors + 1

            errorline = vfh.readline()

    totalassemblybasesinaligns = benchmark_stats["testmattotalcovered"] + benchmark_stats["testpattotalcovered"]
    if totalassemblybasesinaligns > 0:
        phaseerrorrate = totalphasingerrors / totalassemblybasesinaligns
        consensuserrorrate = totalconsensuserrors / totalassemblybasesinaligns
    else:
        phaseerrorrate = 'NA'
        consensuserrorrate = 'NA'

    if consensuserrorrate != "NA":
        if consensuserrorrate > 0:
            consensusqv = int(-10 * math.log(consensuserrorrate, 10))
        else:
            consensusqv = 'Inf'

        if consensuserrorrate + phaseerrorrate > 0:
            qvwithphaseerrors = int(-10 * math.log(consensuserrorrate + phaseerrorrate, 10))
        else:
            qvwithphaseerrors = 'Inf'
    else:
        consensusqv = 'NA'
        qvwithphaseerrors = 'NA'

    generalstatsfile = bedfiles["generalstatsfile"]
    with open(generalstatsfile, "a") as gsfh:
        gsfh.write("\nAssembly accuracy:\n\n")
        gsfh.write("Total errors in alignments of " + args.assembly + " to " + args.benchmark + ": " + str(totalerrors) + "\n")
        gsfh.write("Errors where sequence matches the opposite haplotype: " + str(totalphasingerrors) + "\n")
        gsfh.write("Errors where sequence doesn\'t match the opposite haplotype: " + str(totalconsensuserrors) + "\n")
        gsfh.write("QV without wrong-haplotype errors: " + str(consensusqv) + "\n")
        gsfh.write("QV including wrong-haplotype errors: " + str(qvwithphaseerrors) + "\n")

    return 0

def write_mononuc_stats(mononucstats:dict, bedfiles:dict, benchmark_stats:dict, args):

    totalcoveredmononucs = len(mononucstats.keys())
    totalcorrect = 0
    totalwronglength = 0
    totalphaseswitches = 0
    totalcomplex = 0

    for variantname in mononucstats.keys():
        mononucdict = mononucstats[variantname]
        repeatedbase = mononucdict["base"]
        varianttype = mononucdict["type"]
        reflength = mononucdict["length"]
        newlength = mononucdict["assemblylength"]
        if newlength == -1:
            totalcomplex = totalcomplex + 1
        elif varianttype == "CORRECT":
            totalcorrect = totalcorrect + 1
        elif varianttype == "PHASING":
            totalphaseswitches = totalphaseswitches + 1
        elif varianttype == "CONSENSUS":
            totalwronglength = totalwronglength + 1

    generalstatsfile = bedfiles["generalstatsfile"]
    with open(generalstatsfile, "a") as gsfh:
        gsfh.write("\nMononucleotide run accuracy:\n\n")
        gsfh.write("Total number of true assembly mononucleotide runs of ten or more bases covered by " + args.assembly + " alignments: " + str(totalcoveredmononucs) + "\n")
        perccorrect = round(100*totalcorrect/totalcoveredmononucs, 3)
        gsfh.write("Number of mononucleotide runs correct in assembly: " + str(totalcorrect) + " (" + str(perccorrect) + "%)" + "\n")
        percphaseswitch = round(100*totalphaseswitches/totalcoveredmononucs, 3)
        gsfh.write("Number of mononucleotide runs with assembly alleles matching the alternate haplotype: " + str(totalphaseswitches) + " (" + str(percphaseswitch) + "%)" + "\n")
        percwronglength = round(100*totalwronglength/totalcoveredmononucs, 3)
        gsfh.write("Number of mononucleotide runs with non-alternate alleles of fewer or more of the same base in the assembly: " + str(totalwronglength) + " (" + str(percwronglength) + "%)" + "\n")
        perccomplex = round(100*totalcomplex/totalcoveredmononucs, 3)
        gsfh.write("Number of mononucleotide runs with erroneous alleles other than extensions or contractions: " + str(totalcomplex) + " (" + str(perccomplex) + "%)" + "\n")

    return 0
