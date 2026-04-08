#!/home/tc557/miniconda3/envs/snakemake_new/bin/python
"""
This script makes a multihetsep (mhs) file from a vcf and a series of bed files.
"""

import numpy as np 
import pandas as pd
import argparse
import pdb
import os
import time
import resource

import numpy as np

def collapse_k_clusters_np(seq: np.ndarray, m: int = 2) -> np.ndarray:
    """
    Collapse clusters of heterozygous sites in an integer-encoded sequence.

      mask = -1
      hom  =  0
      het  =  1

    A cluster is: 1 (0{0..m} 1)+
    Only the first 1 in the cluster remains 1.
    All subsequent 1s in the same cluster become 0.

    Any -1 breaks clusters.

    Parameters
    ----------
    seq : np.ndarray, dtype=int
        Array containing {-1, 0, 1}.
    m : int
        Max number of consecutive 0s allowed between 1s for them to count as part of same cluster.

    Returns
    -------
    np.ndarray
        Collapsed sequence (same shape/type as input).
    """
    out = seq.copy()
    in_cluster = False
    t_since_last_k = 0

    for i, val in enumerate(seq):
        if not in_cluster:
            if val == 1:  # start a new cluster
                in_cluster = True
                t_since_last_k = 0
                out[i] = 1
            else:
                # just copy; -1 or 0 cannot start a cluster
                out[i] = val

        else:  # currently inside a cluster
            if val == 1:
                # another K
                if t_since_last_k <= m:
                    # collapse this one → make it hom (0)
                    out[i] = 0
                    t_since_last_k = 0
                else:
                    # too far: start new cluster
                    out[i] = 1
                    t_since_last_k = 0
            elif val == 0:
                # hom, still inside cluster
                out[i] = 0
                t_since_last_k += 1
                if t_since_last_k > m:
                    in_cluster = False
            else:
                # mask; breaks cluster
                out[i] = -1
                in_cluster = False
                t_since_last_k = 0

    return out


def invert_bed(df, chrom_length):
    """
    Given a BED-like DataFrame with columns ['chrom','start','end'],
    return the complement intervals for each chromosome.

    If chrom_sizes is provided (dict: chrom → length),
    the function will include the terminal interval to the chromosome end.
    """
    import pandas as pd
    
    out = []

    # Group by chromosome
    for chrom, sub in df.groupby('chrom'):
        sub = sub.sort_values('start')
        
        prev_end = 0
        
        for _, row in sub.iterrows():
            start, end = row['start'], row['end']
            
            if start > prev_end:
                # gap before this interval
                out.append([chrom, prev_end, start])
            
            prev_end = max(prev_end, end)
        
        # Add gap from last interval to chromosome end
        if prev_end < chrom_length:
            out.append([chrom, prev_end, chrom_length])

    return pd.DataFrame(out, columns=['chrom','start','end'])

def write_mhs(chrom,pos,sspss,gt,filename):
    # pos is index of hets
    # chrom is int
    with open(filename,'w') as f:
        lis=[chrom,pos,sspss,gt]
        for x in zip(*lis):
            f.write("{0}\t{1}\t{2}\t{3}\n".format(*x))
    print(f'\twritten mhs file to {filename}')
    return None

def get_heterozygosity(sequence):
    num_hets = len(sequence[sequence==1])
    num_masks = len(sequence[sequence==-1])
    num_homs = len(sequence[sequence==0])
    print(f'num_hets = {num_hets}')
    print(f'num_homs = {num_homs}')
    print(f'num_masks = {num_masks}')
    est_het = num_hets / (num_hets + num_homs)
    print(f'Est het = {est_het}')
    return None

def test_bed_nonempty(bedfile):
    with open(bedfile,'rb') as f: lines = f.readlines()
    if len(lines) == 0:
        return 0
    else:
        return 1

# Usage : 
# /home/tc557/T2T_PSMC/T2Tsequences/251228_make_mhs.py -iv /home/tc557/rds/hpc-work/T2Tsequences/HG002_all_pap.dip.variantsummary.chr2.0based.txt.gz -negative_beds /home/tc557/rds/hpc-work/T2Tsequences/HG002_all_pap.dip.chr2.indel.bed.gz /home/tc557/rds/hpc-work/T2Tsequences/trd_sdust_test_mappability_on_HG002pri.chr2.sorted.bed -positive_beds /home/tc557/rds/hpc-work/T2Tsequences/HG002_all_pap.dip.bed -chrom 2 -o /tmp/deleteme123

if __name__ == "__main__":
    gstart = time.time()
    parser = argparse.ArgumentParser(description="Set inputs and outputs")
    parser.add_argument('-iv','--in_variants',help='Input variants',required=True,type=str)
    parser.add_argument('-positive_beds','--positive_beds',help='Positive bed file(s). Space delimited if more than one. The positions described are callable, so mask the inverted positions',required=False,type=str, nargs='+',default='')
    parser.add_argument('-negative_beds','--negative_beds',help='Negative bed file(s). Space delimited if more than one. The positions described are not callable, so mask them.',required=True,type=str, nargs='+')
    parser.add_argument('-o','--outfile',help='Output file for mhs',required=True,type=str)
    parser.add_argument('-m','--m',help='Maximum number of 0\'s allowed between 1\'s to be considered in the same cluster. Set -1 for no cleaning',required=False,type=int,default=2)
    parser.add_argument('-chrom','--chrom',help='Chromosome',required=True,type=int)
    parser.add_argument('-truncate','--truncate',help='Truncate start and end by this many base pairs',required=False,type=int,default=50000)
    parser.add_argument('-centromere_start','--centromere_start',help='Start of centromere - mask from here to centromere_end. Set as -1 for no masking',required=False,type=int,default=-1)
    parser.add_argument('-centromere_end','--centromere_end',help='End of centromere - mask from centromere_start to centromere_end. Set as -1 for no masking',required=False,type=int,default=-1)

    args = parser.parse_args()
    zargs = dir(args)
    zargs = [zarg for zarg in zargs if zarg[0]!='_']
    for zarg in zargs:
        print(f'{zarg} is ',end='')
        exec(f'{zarg}=args.{zarg}')
        exec(f'print(args.{zarg})')

    num_negative_masks = len(negative_beds)
    num_positive_masks = len(positive_beds)

    print(f'Number of positive_beds = {num_positive_masks}',flush=True)
    print(f'Number of negative beds = {num_negative_masks}',flush=True)

    print(f'Loading variants...',flush=True)
    variants = pd.read_csv(in_variants,sep=' ',header=None)
    print(f'\tmax length = {variants.iloc[-1,0]}',flush=True)
    possible_lengths = [variants.iloc[-1,0]]
    print(f'Loading negative bed files...',flush=True)
    for bedfile in negative_beds:
        print(f'\t{bedfile}',flush=True)
        try:
            bed = pd.read_csv(bedfile,sep='\t',header=None)
        except:
            print(f'\tWARNING! Failed to load {bedfile} ; empty or wrong format. Skipping...',flush=True)
            continue
        mask = bed[0].str.match(rf"^chr{chrom}($|\D)")
        bed = bed[mask].copy()
        possible_lengths.append(bed.iloc[-1,2])
        print(f'\t\tmax length = {bed.iloc[-1,2]}',flush=True)
    print(f'Loading positive bed files...',flush=True)
    for bedfile in positive_beds:
        print(f'\t{bedfile}',flush=True)
        if test_bed_nonempty(bedfile) == 0:
            print(f'\t\tWARNING! Empty bed file, skipping...',flush=True)
            continue
        bed = pd.read_csv(bedfile,sep='\t',header=None)
        mask = bed[0].str.match(rf"^chr{chrom}($|\D)")
        bed = bed[mask].copy()
        possible_lengths.append(bed.iloc[-1,2])
        print(f'\t\tmax length = {bed.iloc[-1,2]}',flush=True)
    length = max( possible_lengths) + 200
    print(f'Taking max length = {length}',flush=True)
      
    sequence = np.zeros(length,dtype=int)
    genotypes = ['']*length
    
    print(f'Writing variants...')
    for i in range(0,variants.shape[0]):
        sequence[variants.iloc[i,0]] = 1
        genotypes[variants.iloc[i,0]] = variants.iloc[i,1]
    num_hets = len(sequence[sequence==1])
    print(f'\tnum_hets = {num_hets}')
    
    print(f'Writing negative bed files...',flush=True)
    for bedfile in negative_beds:
        num_masks = 0
        print(f'\t{bedfile}',flush=True)
        try:
            bed = pd.read_csv(bedfile,sep='\t',header=None)
        except:
            print(f'\tWARNING! Failed to load {bedfile} ; empty or wrong format. Skipping...',flush=True)
            continue
        bed = pd.read_csv(bedfile,sep='\t',header=None)
        mask = bed[0].str.match(rf"^chr{chrom}($|\D)")
        bed = bed[mask].copy()
        for i in range(0,bed.shape[0]):
            sequence[bed.iloc[i,1]:bed.iloc[i,2]] = -1
            num_masks += bed.iloc[i,2] - bed.iloc[i,1]
        print(f'\t\tnum_masks = {num_masks} = {round(100*num_masks/length,1)}%',flush=True)
    print(f'Writing positive bed files...',flush=True)
    for bedfile in positive_beds:
        num_masks = 0
        print(f'\t{bedfile}',flush=True)
        if test_bed_nonempty(bedfile) == 0:
            print(f'WARNING! Empty bed file, skipping...',flush=True)
            continue
        bed = pd.read_csv(bedfile,sep='\t',header=None)
        mask = bed[0].str.match(rf"^chr{chrom}($|\D)")
        bed = bed[mask].copy()
        num_columns = bed.shape[1]
        bed.columns = ['chrom','start','end'] + ['?']*(num_columns-3)
        invertedbed = invert_bed(bed, length)
        for i in range(0,invertedbed.shape[0]):
            sequence[invertedbed.iloc[i,1]:invertedbed.iloc[i,2]] = -1
            num_masks += invertedbed.iloc[i,2] - invertedbed.iloc[i,1]
        print(f'\t\tnum_masks = {num_masks} = {round(100*num_masks/length,1)}%',flush=True)


    if truncate > 0:
        print(f'Masking first {truncate} and last {truncate} base pairs...',flush=True)
        sequence[0:truncate] = -1
        sequence[-truncate:] = -1
    
    if centromere_start != -1 and centromere_end != -1:
        assert centromere_end > centromere_start, "ERROR! centromere_end must be bigger than centromere_start"
        print(f'Masking centromere; from {centromere_start} - {centromere_end}',flush=True)
        sequence[centromere_start:centromere_end] = -1
        
    num_masks = len(sequence[sequence==-1])
    print(f'num_masks = {num_masks}')
    get_heterozygosity(sequence)

    if m != -1:
        print(f'Cleaning sequence with m = {m}')
        sequence = collapse_k_clusters_np(sequence, m)
        get_heterozygosity(sequence)

    adjacents_ = np.where((sequence[:-1] == 1) & (sequence[1:] == 1))[0]
    if m != -1:
        assert len(adjacents_) == 0, "ERROR! There should be no adjacent hets, but apparently there are. This suggests an error."
    else:
        if len(adjacents_) != 0:
            print(f'WARNING! There are adjacent {len(adjacents_)} hets (setting m>-1 should remove this). Not aborting.',flush=True)

    het_positions = np.where(sequence==1)[0]
    gts = [genotypes[i] for i in het_positions]
    sspss = [] # (number of callable ) sites since previous segregating sites
    prev = 0
    for index,i in enumerate(het_positions):
        # if (index%1000)==0: print(f'on het {index} out of {num_hets}',flush=True)
        zsspss = int(i - prev + sequence[prev:i].sum() + 1)
        sspss.append( zsspss )
        prev = i + 1
    sspss = np.array(sspss)
    assert sspss.min() >= 1, f"Min sspss can not be less than one, but it is {sspss.min()}. There must be an error. Aborting."
    chrom_str = [f'chr{str(chrom)}']*len(sspss)
    # gts =  [string_gts[(jj%length_stringgts):((jj+2)%length_stringgts[10:])] for jj in range(len(het_positions))]
    write_mhs(chrom_str,het_positions,sspss,gts,outfile)

    usage = resource.getrusage(resource.RUSAGE_SELF)
    max_memory_usage = usage.ru_maxrss  # memory used in bytes
    gdone = time.time()
    print(f'Done!',flush=True)
    print(f'\tTime taken: {round(gdone - gstart,2)}s',flush=True)
    print(f'\tMax memory used: {max_memory_usage / 1024:.2f}MB', flush=True)
    
