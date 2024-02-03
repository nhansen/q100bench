import re
import pysam
import pybedtools
from collections import namedtuple
from q100bench import seqparse

# create namedtuple for bed intervals:
bedinterval = namedtuple('bedinterval', ['chrom', 'start', 'end', 'name', 'rest']) 


