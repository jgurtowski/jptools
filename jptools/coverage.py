from __future__ import print_function

import sys
from itertools import groupby, imap, chain, izip, repeat, count
from functools import partial
from operator import attrgetter, itemgetter


from jbio.io.file import iterator_over_file
from jbio.io.blast import record_iterator as blast_record_iterator
from jbio.coverage import coverage_array_from_ranges
from jbio.coverage import get_marked_ranges
from jbio.functional import compose
from jbio.alignment import best_scoring_non_overlapping

def coverage_from_blast6():

    if not len(sys.argv) == 2:
        sys.exit("coverage_from_blast6 in.blast6")

    raw_alignment_it = blast_record_iterator(iterator_over_file(sys.argv[1]))
    lno = partial(best_scoring_non_overlapping,
                  attrgetter("qstart"), 
                  attrgetter("qend"),
                  attrgetter("bitscore"))

    q_filt_alignment_it = chain.from_iterable(
        imap(compose(lno, itemgetter(1)), 
             groupby(raw_alignment_it, 
                     attrgetter("qname"))))
    
    #read all alignments into memory
    ref_sorted_alignments = sorted(q_filt_alignment_it, 
                                   key=attrgetter("sname"))
    
    for reference,alignments in groupby(ref_sorted_alignments,
                                        attrgetter("sname")):
        alignments = list(alignments)
        ref_len = alignments[0].slen
        
        blast_start_getter = lambda a: a.sstart-1
        blast_end_getter = lambda a: a.send-1
        cov_arr = coverage_array_from_ranges(alignments, ref_len,
                                             blast_start_getter,
                                             blast_end_getter)
        filter(print, izip(count(1),cov_arr))
        #mark the regions with 0 coverage
        zerocov = map(lambda x: 1 if x==0 else 0, cov_arr)
        zerocov_regions = get_marked_ranges(zerocov)
        
        region_printer = compose(print,lambda (x,(y,z)) : "\t".join(map(str,[x,y,z])))
        
        #filter(region_printer, izip(repeat(reference), zerocov_regions))
        
        ##Get Low ID regions
        ranges_w_id = imap(compose(lambda (x,y,i) : (x-1,y-1,i),
                                   attrgetter("sstart","send","pctid")), alignments)

        pct_arr = coverage_array_from_ranges(alignments, ref_len,
                                             blast_start_getter,
                                             blast_end_getter,
                                             lambda r, (o_pid,o_cnt): (r.pctid+o_pid, o_cnt+1), 
                                             (0,0))
        lowid = map(lambda (c_pid,cnt): 1 if cnt != 0 and c_pid/cnt < 95.0 else 0,
                    pct_arr)
        
        lowid_regions = get_marked_ranges(lowid)
        
        #filter(region_printer, izip(repeat(reference), lowid_regions))
