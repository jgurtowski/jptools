from __future__ import print_function
from itertools import chain, imap, tee, izip
from functools import partial
from operator import attrgetter, itemgetter

from jbio.io.file import iterator_over_file, record_to_string
from jbio.io.blast import record_iterator as blast_record_iterator
from jbio.functional import compose
from jbio.alignment import best_scoring_non_overlapping
from jbio.alignment import group as alignment_grouper, LIS

def blast6filter_main(cmdline_args = None):
    
    if not cmdline_args:
        import sys
        cmdline_args = sys.argv

    if not len(cmdline_args) == 3:
        return "blast6filter q/r input.blast6 -- Make sure input is sorted by q/r first"
    
    task,infile = cmdline_args[1:3]

    fileit = iterator_over_file(infile)

    alignment_getter = blast_record_iterator(fileit)
    
    grouped_alns = alignment_grouper(attrgetter("sname") if task == "r" else attrgetter("qname")
                                     ,alignment_getter)
    if task == "r":
        best = imap(compose(partial(map,itemgetter(2)),partial(LIS,attrgetter("sstart"), attrgetter("send"), None)),
                    grouped_alns)
    else:
        best = imap(partial(best_scoring_non_overlapping,attrgetter("qstart"),
                            attrgetter("qend"),
                            attrgetter("bitscore")),
                    grouped_alns)
    
    filter(print,imap(record_to_string, 
                      chain.from_iterable(best)))
    
    
