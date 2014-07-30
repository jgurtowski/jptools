from __future__ import print_function

import sys
from itertools import groupby, imap, chain
from functools import partial
from operator import attrgetter, itemgetter

from jbio.io.file import iterator_over_file, line_record_iterator
from jbio.io.blast import Blast6Record, Blast6Types
from jbio.coverage import coverage_array_from_ranges
from jbio.coverage import get_marked_ranges
from jbio.functional import compose
from jbio.alignment import best_scoring_non_overlapping

def coverage_from_blast6():

    if not len(sys.argv) == 2:
        sys.exit("coverage_from_blast6 in.blast6")

    raw_alignment_it = line_record_iterator(Blast6Record, Blast6Types, 
                                            iterator_over_file(sys.argv[1]))
    lno = partial(best_scoring_non_overlapping,
                  attrgetter("qstart"), 
                  attrgetter("qend"),
                  attrgetter("bitscore"))

    q_filt_alignment_it = chain.from_iterable(
        imap(compose(lno, itemgetter(1)), 
             groupby(raw_alignment_it, 
                     attrgetter("qname"))))
    
    ref_sorted_alignments = sorted(q_filt_alignment_it, 
                                   key=attrgetter("sname"))
    
    for reference,alignments in groupby(ref_sorted_alignments,
                                        attrgetter("sname")):
        alignments = list(alignments)
        ref_len = alignments[0].slen
        ranges = imap(attrgetter("sstart","send"), alignments)
        cov_arr = coverage_array_from_ranges(ranges, ref_len)
        #mark the regions with 0 coverage
        zerocov_regions = map(lambda x: 1 if x==0 else 0, cov_arr)
        print(get_marked_ranges(zerocov_regions))

        
