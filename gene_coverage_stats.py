#!/usr/bin/env python
from __future__ import print_function

##
#Determines how many genes are represented intact
#within an assembly
#
import sys
from itertools import groupby, imap
from operator import attrgetter, add
import functools

from jbio.io.blast import record_iterator as blast_record_iterator
from jbio.io.file import iterator_over_file
from jbio.alignment import alignment_functions

if not len(sys.argv) == 2:
    sys.exit("gene_coverage_stats.py input.blast6.q")

infile = sys.argv[1]

blast_records = blast_record_iterator(iterator_over_file(infile))

afuncs = alignment_functions(attrgetter("qstart"), attrgetter("qend"))

for qname, hsps in groupby(blast_records, attrgetter("qname")):
    hsps = list(hsps)
    first_hsp = hsps[0]
    
    hsp_lens = map(afuncs.len_aln, hsps)
    hsp_pct_of_query = map(lambda l : l / float(first_hsp.qlen) * 100.0, hsp_lens)
    total_query_coverage = sum(hsp_lens) / float(first_hsp.qlen) * 100.0
    hsp_contig_names = map(attrgetter("sname"), hsps)
    number_different_contigs = len(set(sorted(hsp_contig_names)))
    
    p_items = [qname, total_query_coverage, number_different_contigs, zip(hsp_pct_of_query, hsp_contig_names)]
    print("\t".join(imap(str, p_items)))
    


