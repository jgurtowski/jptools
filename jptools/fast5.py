import os
import h5py
import logging
from operator import itemgetter

from itertools import imap, izip, ifilter

def fast5ToFasta(f5files, outfh, logfh):
    '''Takes fast5files and outputs to output_file and log_file'''
    opj = os.path.join
    
    def fqstrToArr(fqstr):
        arr = fqstr.strip().split()
        if not len(arr) == 4:
            logfh("Fastq string looks corrupt %s \n" % arr[0])
            return None
        tname = arr[0][1:]
        arr[0] = tname
        tdesc = arr[2][1:]
        arr[2] = tdesc
        return arr

    for f5file in f5files:
        with h5py.File(f5file) as f5fh:
            logfh.write("Fast5File %s\n" % f5file)
            
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
        
            fq_datasets =  map(lambda x : f5fh.get(x), fq_grp_names)
        
            fq_arrs = imap(lambda x : fqstrToArr(x[()]), fq_datasets)
        
             #(bcg_suffix, [fastq_arr])
            filt_arrs = ifilter(lambda x : not itemgetter(1)(x) == None,  izip(bcg_suffix, fq_arrs))
            fastas = imap(lambda x : (x[1][0]+"_"+exp_start+"_"+x[0], x[1][1]), filt_arrs)
                    
            for fasta in fastas:
                outfh.write( "\n".join([">"+fasta[0], fasta[1]]) )
                outfh.write( "\n")
