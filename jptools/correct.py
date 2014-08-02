#!/usr/bin/env python

import sys

from operator import attrgetter
from itertools import groupby, imap
from collections import namedtuple

from jbio.io.blast import Blast6SeqRecord, Blast6SeqTypes

from jbio.misc import reverse_complement
from jbio.sequence import fasta_iterator
from jbio.io.file import iterator_over_file, line_record_iterator
from jbio.io.blast import Blast6Record, Blast6Types
from jbio.functional import compose

from jbio.log import logger

from pbtools.pbdagcon.c_aligngraph import AlnGraph, convert_mismatches
#from pbtools.pbdagcon.aligngraph import AlnGraph, convert_mismatches
from pbtools.pbdagcon.q_sense import output_dag_info


TOO_MANY_ALIGNMENTS = 3000

def correct_oxford(reads_fn=None, alignments_fn=None):
    '''Corrects oxford reads'''
    
    log = logger(sys.stderr)
    
    if not reads_fn or not alignments_blast6_fn:
        if not len(sys.argv) == 3:
            sys.exit("correct.py raw_reads.fa alignments.blast6")
        (reads_fn,alignments_fn) = sys.argv[1:3]

        log("Reading raw reads into memory")
        #just put all reads in memory
        fastas = compose(fasta_iterator, iterator_over_file)(reads_fn)
        raw_reads = dict(map(attrgetter("name","seq"), fastas))

        log("Reading raw reads DONE :)")

        #The alignments need to be sorted by the long read name (second column)
        alignment_it = line_record_iterator(Blast6SeqRecord, Blast6SeqTypes,
                                            iterator_over_file(alignments_fn))
        
        important_field_getter = attrgetter("qname","sname","qstart","qend",
                                            "sstart","send", "qseq", "sseq")
                                            
        for readname, alignments in groupby(alignment_it, attrgetter("sname")):
            log("Working on %s" % readname)
            
            raw_read_seq = raw_reads.get(readname)
            if not raw_read_seq:
                log("Can not find sequence for %s" % readname)
                continue

            log("Raw Read Length: %d" % len(raw_read_seq))    
            g = AlnGraph(raw_read_seq)

            alignments = imap(important_field_getter, alignments)
            num_alignments = 0
            for qname,sname,qstart,qend,sstart,send,qseq,sseq in alignments:

                #blast alignments are one based, convert to 0 based
                (qstart, qend) = (qstart-1, qend-1)
                (sstart, send) = (sstart-1, send-1)

                #reverse complement, must switch the alignment strings
                if send < sstart:
                    (qseq, sseq) = tuple(map(reverse_complement, [qseq,sseq]))
                    send, sstart = sstart,send
                    
                (qseq, sseq) = convert_mismatches(qseq,sseq)
                try:
                    alignment_tuple =((qstart, qend, qseq),
                                      (sstart, send, sseq), qname) 
                    g.add_alignment( alignment_tuple)
                except Exception as e:
                    log("Add Alignmented Error: %s" % e)
                    continue
                if num_alignments > TOO_MANY_ALIGNMENTS:
                    break
                
                num_alignments += 1

            log("Processed Alignments: %d" % num_alignments)
            if num_alignments > TOO_MANY_ALIGNMENTS:
                log("Too Many Alignments, Skipping")
                continue
            
            log("Generating Consensus")
            consensus = g.generate_all_consensus(min_cov=0)[0]
            log("Consensus Length %d" % len(consensus[0]))
            log("%s Done\n\n" % readname)

            log("Output dag info")
            output_dag_info(g, "g.info")

            print ">"+readname+"_consensus"
            print consensus[0]


