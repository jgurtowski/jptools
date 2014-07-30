from __future__ import print_function
from itertools import chain, imap, tee, izip
from functools import partial
from operator import attrgetter

from jbio.io.file import iterator_over_file, record_to_string
from jbio.io.blast import record_iterator as blast_record_iterator
from jbio.alignment import best_scoring_non_overlapping
from jbio.alignment import group as alignment_grouper

def blast6qfilter_main(cmdline_args = None):
    
    if not cmdline_args:
        import sys
        cmdline_args = sys.argv

    if not len(cmdline_args) == 2:
        return "blast6filter input.blast6"
    
    infile = cmdline_args[1]

    fileit = iterator_over_file(infile)

    alignment_getter = blast_record_iterator(fileit)
    
    q_grouped_alns = alignment_grouper(attrgetter("qname"),
                                       alignment_getter)
    
    
    q_grouped_best = imap(partial(best_scoring_non_overlapping,attrgetter("qstart"),
                                                  attrgetter("qend"),
                                                  attrgetter("bitscore")),
                          q_grouped_alns)
    
    filter(print,imap(record_to_string, 
                      chain.from_iterable(q_grouped_best)))
    
