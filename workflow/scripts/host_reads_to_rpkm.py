#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys


sys.stderr = open(snakemake.log[0], "w")


import pandas as pd


def main(cov, mapping_stats, map):
    # get total reads
    with open(mapping_stats, 'r') as fi:
        for line in fi:
            l = line.split("\t")
            if l[0] == "SN" and l[1] == "raw total sequences:":
                miReads = int(l[2])/1000000
    
    # Adding organism info to coverage table, then collapse and finally calculate rpkm
    map = pd.read_csv(map, sep=" ", names=['id', 'genus', 'species'])
    cov = pd.read_csv(cov, sep="\t"
        ).merge(
            map, left_on="#rname", right_on='id', how='left'
        ).groupby(
            ['genus', 'species']
        ).sum().reset_index()
    cov['RPM'] = cov.apply(lambda x: x['numreads']/miReads, axis = 1)
    cov['RPKM'] = cov.apply(lambda x: x['RPM']/(x['endpos']/1000), axis = 1)
    cov = cov.sort_values('RPKM', ascending=False)
    cov['organism'] = cov['genus'].map(str) + ' ' + cov['species'].map(str)
    cov[['organism', 'RPKM']].to_csv(rpkm, sep="\t", index=False)


if __name__ == '__main__':
    main(snakemake.input['cov'],
         snakemake.input['mapping_stats'],
         snakemake.input['map'],
         snakemake.output['rpkm'])
