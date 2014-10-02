#!/usr/bin/env python

import sys

from operator import itemgetter, attrgetter
from itertools import imap, ifilter, izip
from functools import partial
from urllib import unquote

from jbio.io.file import iterator_over_file
from jbio.gff import record_iterator as gff_iterator
from jbio.fasta import record_iterator as fasta_iterator

if not len(sys.argv) == 3:
    print "gene_fasta.py input.fa input.gff"
    sys.exit(1)

#FIELDS = ["ID","Alias","orf_classification","gene","Note"]
FIELDS = ["ID","Note"]

fa_fn,gff_fn  = sys.argv[1:3]

#read fasta records into memory
def fasta_clean_getter(fasta_entry):
    name = fasta_entry.name.split()[0]
    return (name, fasta_entry.seq)

fasta_records = dict(imap(fasta_clean_getter,fasta_iterator(iterator_over_file(fa_fn))))

gene_entries = ifilter(lambda x: x.feature == "gene", 
                       gff_iterator(iterator_over_file(gff_fn)))

for gene_record in gene_entries:
    attrs = dict(map(lambda x: x.split("="), gene_record.attribute.split(";")))
    header = ">" + attrs["Name"]
    fields = FIELDS
    field_getter_func = lambda x : unquote(attrs.get(x,"None")) if x =="Note" else attrs.get(x,"None")
    field_getter = imap(field_getter_func, fields)
    header += " " + " ".join(imap(lambda fv: "[%s=%s]" % fv, izip(fields, field_getter)))
    
    start, end = gene_record.start-1, gene_record.end-1
    seq = fasta_records[gene_record.seqname][start:end+1]
    print header
    print seq



    
