from __future__ import print_function

import sys

from itertools import chain, imap, tee, izip, cycle, ifilter
from functools import partial
from operator import attrgetter, itemgetter

from jbio.io.file import iterator_over_file, record_to_string
from jbio.io.blast import record_iterator as blast_record_iterator
from jbio.functional import compose, zipmap
from jbio.alignment import group as alignment_grouper, alignment_functions
from jbio.log import logger


log = logger(sys.stderr)

def blast6filter_main(cmdline_args = None):
    
    if not cmdline_args:
        import sys
        cmdline_args = sys.argv

    if not len(cmdline_args) == 3:
        return "blast6filter q/r_cons/r_noover input.blast6 -- Make sure input is sorted by q/r first"
    
    task,infile = cmdline_args[1:3]

    fileit = iterator_over_file(infile)

    alignment_getter = blast_record_iterator(fileit)

    #
    #grouped_alns = alignment_grouper(attrgetter("sname"), alignment_getter)
    #aln_funcs = alignment_functions(attrgetter("sstart"), attrgetter("send"))
    #filter(compose(print,record_to_string), chain.from_iterable(imap(partial(aln_funcs.greedy_repeat_filter,final_sort_key=attrgetter("pctid")), grouped_alns)))
    #
    #sys.exit(1)
    
    if task.startswith("r"):
        grouped_alns = alignment_grouper(attrgetter("sname"), alignment_getter)
        aln_funcs = alignment_functions(attrgetter("sstart"), attrgetter("send"))
        score_func = aln_funcs.score_getter_matching_consensus_estimated
        greedy_repeat_filt = partial(aln_funcs.greedy_repeat_filter, final_sort_key=attrgetter("pctid"))
        def remove_self(alns):
            a = list(alns)
            log("Remove Self: Working on %d" % len(alns))
            filtered = filter(lambda y: not y.qname == y.sname, a)
            log("Remove Self: Filtered alignments: %d" % len(filtered))
            return filtered
            
        lis = compose(partial(aln_funcs.LIS, score_func), aln_funcs.remove_contained, greedy_repeat_filt, remove_self)
        best = imap(lis, grouped_alns)
        if task == "r_noover":
            score_func = aln_funcs.score_getter_penalize_overlap_estimated
            lis = compose(partial(aln_funcs.LIS, score_func), aln_funcs.remove_contained)
            best = imap(compose(partial(map,itemgetter(2)),lis), grouped_alns)        
        if task =="r_experimental":
            lis = compose(greedy_repeat_filt, remove_self)
            best = imap(lis, grouped_alns)

    else:
        grouped_alns = alignment_grouper(attrgetter("qname"), alignment_getter)
        aln_funcs = alignment_functions(attrgetter("qstart"), attrgetter("qend"))
        lis = compose(partial(aln_funcs.LIS,aln_funcs.score_getter_penalize_overlap_estimated), aln_funcs.remove_contained)
        best = imap(compose(partial(map,itemgetter(2)),lis), grouped_alns)        
        
    

    
    filter(print,imap(record_to_string, 
                      chain.from_iterable(best)))
    
    
