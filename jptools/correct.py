#!/usr/bin/env python

import sys

from operator import attrgetter
from itertools import groupby, imap
from collections import namedtuple

from jbio.io.blast import Blast6SeqRecord, Blast6SeqTypes

from ectools.misc import reverse_complement
from ectools.seqio import fastaIterator
from ectools.nucio import fileIterator, lineRecordIterator
from ectools.cov import getCoverageFromAlignments
from ectools.log import logger

from pbtools.pbdagcon.c_aligngraph import AlnGraph, convert_mismatches
from pbtools.pbdagcon.q_sense import output_dag_info

log = logger(sys.stderr)

def correct_oxford(reads_fn=None, alignments_fn=None):
    '''Corrects oxford reads'''
    
    if not reads_fn or not alignments_blast6_fn:
        if not len(sys.argv) == 3:
            sys.exit("correct.py raw_reads.fa alignments.blast6")
        (reads_fn,alignments_fn) = sys.argv[1:3]

        log("Reading raw reads into memory")
        #just put all reads in memory
        raw_reads = dict(map(attrgetter("name","seq"), 
                             fileIterator(reads_fn, fastaIterator)))

        log("Reading raw reads DONE")

         #The alignments need to be sorted by the long read name (second column)
        record_it = lambda fh : lineRecordIterator(fh, Blast6SeqRecord, Blast6SeqTypes)
        file_it = fileIterator(alignments_fn, record_it)

        important_field_getter = attrgetter("qname","sname","qstart","qend",
                                            "sstart","send", "qseq", "sseq")

        for readname, alignments in groupby(file_it, attrgetter("sname")):
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
                    g.add_alignment( ((qstart, qend, qseq),
                                      (sstart, send, sseq), qname))
                except Exception as e:
                    log("Add Alignmented Error: %s" % e)
                    continue

                num_alignments += 1

            log("Processed Alignments: %d" % num_alignments)
            log("Generating Consensus")
            consensus = g.generate_consensus(min_cov=0)[0]
            log("Consensus Length %d" % len(consensus))
            log("%s Done\n\n" % readname)

            print ">"+readname+"_consensus"
            print consensus


