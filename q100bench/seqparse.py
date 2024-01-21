import re
import pysam
import pybedtools
from collections import namedtuple
import time

# create namedtuple for bed intervals:
bedinterval = namedtuple('bedinterval', ['chrom', 'start', 'end', 'name', 'rest']) 

def write_genome_bedfile(queryobj, args, genome_bedfile=None )->pybedtools.BedTool:
    genomebedstring = ""
    with open(genome_bedfile, "w") as gfh:
        for scaffold in queryobj.references:
            scaffoldlength = queryobj.get_reference_length(scaffold)
            scaffstring = scaffold + "\t0\t" + str(scaffoldlength - 1) + "\n"
            gfh.write(scaffstring)
            genomebedstring += scaffstring

    return pybedtools.BedTool(genomebedstring, from_string = True)

def find_all_ns(queryobj, args, n_bedfile=None, atgc_bedfile=None )->list:

    # did the command-line arguments specify a pre-existing bedfile of N-stretch locations?
    user_n_file = args.n_bedfile
    if n_bedfile:
        nfh = open(n_bedfile, "w")
    if atgc_bedfile:
        atgcfh = open(atgc_bedfile, "w")

    if not user_n_file:
        p = re.compile("N+")
    else:
        ofh = open(user_n_file, "r")

    n_intervals = []
    contig_intervals = []
    n_interval_dict = {}
    # if an n interval file has been specified in the command line options, read it:
    if user_n_file:
        nbed_line = ofh.readline()
        while nbed_line:
            [chrom, start, end, gapname] = nbed_line.split("\t")
            n_intervals.append(bedinterval(chrom=chrom, start=int(start), end=int(end), name=gapname, rest=""))
            if chrom in n_interval_dict.keys():
                n_interval_dict[chrom].append(bedinterval(chrom=chrom, start=int(start), end=int(end), name=gapname, rest=""))
            else:
                n_interval_dict[chrom] = [bedinterval(chrom=chrom, start=int(start), end=int(end), name=gapname, rest="")]
            nbed_line = ofh.readline()
        ofh.close()
        for ref in queryobj.references:
            if ref not in n_interval_dict.keys():
                n_interval_dict[ref] = []
            contignum = 1
            contigstart = 0
            for interval in n_interval_dict[ref]:
                interval_name = ref + "." + str(contignum)
                contig_intervals.append(bedinterval(chrom=ref, start=contigstart, end=interval.start, name=interval_name, rest=""))
                contignum = contignum + 1
                contigstart = interval.end
            refend = queryobj.get_reference_length(ref)
            interval_name = ref + "." + str(contignum)
            contig_intervals.append(bedinterval(chrom=ref, start=contigstart, end=refend, name=interval_name, rest=""))
    else:
        for ref in queryobj.references:
            contignum = 1
            contigstart = 0
            findstring = 'N' * args.minns
            if n_bedfile:
                #startfindn = time.time()
                chromseq = queryobj.fetch(ref).upper()
                refend = queryobj.get_reference_length(ref)
                start = chromseq.find(findstring)
                while start != -1:
                    end = start + args.minns
                    while end < refend and chromseq[end] == 'N':
                        end = end + 1
                    if end - start >= args.minns:
                        if n_bedfile:
                            gapname = "N." + ref + "." + str(contignum)
                            nfh.write(ref + "\t" + str(start) + "\t" + str(end) + "\t" + gapname + "\n")
                            n_intervals.append(bedinterval(chrom=ref, start=start, end=end, name=gapname, rest=""))
                        contigend = start
                        contigname = ref + "." + str(contignum)
                        contig_intervals.append(bedinterval(chrom=ref, start=contigstart, end=contigend, name=contigname, rest=""))
                        contignum = contignum + 1
                        contigstart = end
                    start = chromseq.find(findstring, end, refend)

                #endfindn = time.time()
                #print("Find ns time new way for chrom %s: %s" % (ref, endfindn - startfindn))

                contigname = ref + "." + str(contignum)
                contig_intervals.append(bedinterval(chrom=ref, start=contigstart, end=refend, name=contigname, rest=""))
        if n_bedfile:
            nfh.close()
    if atgc_bedfile:
        for interval in contig_intervals:
            atgcfh.write(interval.chrom + "\t" + str(interval.start) + "\t" + str(interval.end) + "\t" + interval.name + "\n")
        atgcfh.close()

    return [contig_intervals, n_intervals]

def revcomp(seq:str) -> str:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'M': 'N', 'K': 'N', 'R': 'N', 'W': 'N', 'Y': 'N'}
    bases = list(seq)
    bases = bases[::-1]
    bases = [complement[base] for base in bases]

    return ''.join(bases)

