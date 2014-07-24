
import sys
from itertools import groupby, imap
from operator import attrgetter


from ectools.m4io import getRawAlignments, getAlignments
from ectools.blastio import Blast6Record, Blast6Types
from jbio.coverage import coverage_array_from_ranges
from jbio.coverage import get_marked_ranges

def coverage_from_blast6():

    if not len(sys.argv) == 2:
        sys.exit("coverage_from_blast6 in.blast6")
        
    alignment_it = getRawAlignments(sys.argv[1], Blast6Record, Blast6Types)
    sorted_alignments = sorted(alignment_it, key=attrgetter("sname"))
    
    
    for reference,alignments in groupby(sorted_alignments,
                                        attrgetter("sname")):
        alignments = list(alignments)
        ref_len = alignments[0].slen
        ranges = imap(attrgetter("sstart","send"), alignments)
        cov_arr = coverage_array_from_ranges(ranges, ref_len)
        print "hi"
        #mark the regions with 0 coverage
        zerocov_regions = map(lambda x: 1 if x==0 else 0, cov_arr)
        print get_marked_ranges(zerocov_regions)


