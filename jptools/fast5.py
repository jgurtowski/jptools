import os
import h5py
import logging
import sys
from operator import itemgetter
from functools import partial
from itertools import imap, izip, ifilter, chain

from jbio.fastq import record_iterator as fastq_iterator

def fast5ToFasta_main():
    if not len(sys.argv) >= 2:
        sys.exit("fast5ToFasta [directory] [in.fast5 [in2.fast5]]")
    
    in_files = sys.argv[1:]
    if os.path.isdir(in_files[0]):
        path = in_files[0]
        fds = os.listdir(path)
        in_files = map(partial(os.path.join, path), ifilter(lambda x: x.endswith(".fast5"), fds))
    
    fast5ToFasta(in_files, sys.stdout, sys.stderr)

def fast5ToFasta(f5files, outfh, logfh):
    '''Takes fast5files and outputs to output_file and log_file'''
    opj = os.path.join

    for f5file in f5files:
        logfh.write("Fast5File %s\n" % f5file)
        try:
            with h5py.File(f5file,'r') as f5fh:
                #get experiment start time
                tid_grp = f5fh['/UniqueGlobalKey/tracking_id']
                exp_start = tid_grp.attrs.get("exp_start_time")
        
                #get the different fastq entries (could be template/complement/2D all or some)
                basecall_root = '/Analyses/Basecall_2D_000'
                base_call = f5fh[basecall_root]
        
                basecalled_grps = filter( lambda x: x.startswith("BaseCalled_"), base_call.keys())
            
                bcg_suffix = map( lambda x : x.split("_")[1], basecalled_grps)
                logfh.write("Found %d BaseCall Groups : %s\n" % (len(bcg_suffix), ",".join(bcg_suffix)))
        
                fq_grp_names = map( lambda x : opj(opj(basecall_root, x), "Fastq"), basecalled_grps)

                fq_datasets =  map(lambda y : f5fh.get(y), fq_grp_names)

                fastq_arrays = imap(lambda d: d[()].strip().split("\n"), fq_datasets)
                
                fastq_objs = fastq_iterator( chain.from_iterable(fastq_arrays))
                
                for fastq_obj in fastq_objs:
                    outfh.write( ">" + fastq_obj.name + "\n")
                    outfh.write( fastq_obj.seq + "\n")
        except Exception as e:
            logfh.write("Error : %s" % str(e))
            pass
