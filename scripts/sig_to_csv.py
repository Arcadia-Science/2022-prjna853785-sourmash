#!/usr/bin/env python3
"""
Given a sourmash signature, extract the minhash integers and abund info
to a csv file with the basename of the sig file name as the name of the 
abund column.
python sig_to_csv.py 31 signature csv
"""

from sourmash import signature
import pandas as pd
import os
import sys
import argparse

def main():
    p = argparse.ArgumentParser()
    p.add_argument('ksize')         # kmer size to use; this is brittle as it assumes each signature has sketches with k = 21, k = 31, k = 51 that occur in that order. This is true for this workflow, but should be coded in a more robust way.
    p.add_argument('signature')       # sourmash signature
    p.add_argument('output')          # output csv file name
    args = p.parse_args()

    # load the signature from disk
    sigfp = open(args.signature, 'rt')
    siglist = list(signature.load_signatures(sigfp))
    if args.ksize == "21":
        loaded_sig = siglist[0]
    elif args.ksize == "31":
        loaded_sig = siglist[1]
    elif args.ksize == "51":
        loaded_sig = siglist[2]
    else:
        print("unexpected kmer size, exiting")
    
    mins = loaded_sig.minhash.hashes.keys() # get minhashes
    name = loaded_sig.name # use the signature name as the column name
    df = pd.DataFrame(mins, columns=[name]) # set the column name
    df.to_csv(args.output, index = False) # write to a csv

if __name__ == '__main__':
    sys.exit(main())
