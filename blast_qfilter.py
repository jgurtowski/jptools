#!/usr/bin/env python

import sys

from Bio.Blast.NCBIXML import parse as parseBlast
from ectools.misc import reverse_complement
from ectools.seqio import fastaIterator
from ectools.cov import getCoverageFromAlignments

from pbtools.pbdagcon.c_aligngraph import AlnGraph, convert_mismatches
from pbtools.pbdagcon.q_sense import output_dag_info

if not len(sys.argv) == 3:
    sys.exit("blast_qfilter.py read.fa input.xml")

printer = False

with open(sys.argv[1]) as fafh:
    ref_read = list(fastaIterator(fafh))[0]

g = AlnGraph(ref_read.seq)

read_length = len(ref_read.seq)
alignment_positions = []

def remove_ref_indels(qAln, sAln):
    '''taken from c_utils in pbdagcon'''
    new_aligns = []
    for b1, b2 in zip(qAln,sAln):
        if b1 == "-":
            b1 = b2
        if not b2 == "-":
            new_aligns.append( (b1, b2) )
    (q,s) = zip(*new_aligns)

    return ("".join(q), "".join(s))

with open(sys.argv[2]) as fh:
    for record in parseBlast(fh):
        #print record.query
        for alignment in record.alignments:
            #print alignment.hit_id, alignment.hit_def
            for hsp in alignment.hsps:
                if hsp.expect > 1e-9:
                    continue
                (qAln, sAln) = (str(hsp.query), str(hsp.sbjct))
                (sbjS, sbjE) = (hsp.sbjct_start-1, hsp.sbjct_end-1)
                (qS, qE) = (hsp.query_start-1, hsp.query_end -1)
                if (sbjE < sbjS):
                    (qAln, sAln) = (reverse_complement(qAln), reverse_complement(sAln))
                    t = sbjS
                    sbjS = sbjE
                    sbjE = t

                #(qAln,sAln) = remove_ref_indels(qAln, sAln)
                
                alignment_positions.append( (sbjS, sbjE) )
                if "28313:13381" in str(record.query):
                    sys.stderr.write("%s,%s,%s \n" %(qS, qE, qAln))
                    sys.stderr.write("%s,%s,%s \n" %(sbjS,sbjE, sAln))
                    (qAln, sAln) = convert_mismatches(qAln, sAln)
                    sys.stderr.write("\n\n")
                    sys.stderr.write("%s,%s,%s \n" %(qS, qE, qAln))
                    sys.stderr.write("%s,%s,%s \n" %(sbjS,sbjE, sAln))
                else:
                    (qAln, sAln) = convert_mismatches(qAln, sAln)
                g.add_alignment( ((qS, qE, qAln), 
                                  (sbjS, sbjE, sAln)), record.query)

                if printer:
                    print hsp.align_length, hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, hsp.strand
                    print hsp.query
                    print hsp.match
                    print hsp.sbjct

#print getCoverageFromAlignments(alignment_positions, read_length)
print ">consensus_read"
print g.generate_consensus(min_cov=0)[0]

output_dag_info(g, "dag.info")
