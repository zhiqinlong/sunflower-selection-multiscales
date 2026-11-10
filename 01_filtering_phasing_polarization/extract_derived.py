#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 22:14:18 2026

@author: zhiqinlong
"""
#p_anc: the probability of major allele being ancestral
import pandas as pd
import argparse
parser= argparse.ArgumentParser(
	prog = "extract_derived.py",
	description = "make a file labbled which allele is ancestral")
parser.add_argument('-com',"--outgroup_combination",action='store',help='the path of different combination of outgroup specis')
args = parser.parse_args()

pro_major_df=pd.read_csv(f"{args.outgroup_combination}.est_sfs.output.p_anc.txt",sep=r'\s+',header=None,skiprows=8,index_col=False )
pro_major=pro_major_df.iloc[:,2].tolist()
major_allele_df= pd.read_csv(f"{args.outgroup_combination}_major.txt", sep='\t')
major_allele_df['ID'] = major_allele_df['CHROM'] + ':' + major_allele_df['POS'].astype(str)
orders=['A', 'C', 'G', 'T']
f_out=open(f"{args.outgroup_combination}.derived_out.txt",'w') 
f_out.write('\t'.join(major_allele_df.columns.tolist())+ "\t"+"ancestral"+"\n")
count = sum(x > 0.5 for x in pro_major)
print(count)  
options_an = ['REF', 'ALT']
for idx, row in major_allele_df.iterrows():
    ancestral=''
    if row['major']=='nomajor':
        if orders.index(row['REF'])<orders.index(row['ALT']):
            major_a='REF'
            if pro_major[idx]>0.5:        
                ancestral=major_a
            else:      
                ancestral='ALT'
        else:
            if pro_major[idx]>0.5:      
                major_a='ALT'
                ancestral=major_a
            else:      
                ancestral='REF'
    else:
        if pro_major[idx]>0.5:
            ancestral=row['major'] 
        else:
            other = [a for a in options_an if a != row['major']][0]
            ancestral = other 
    line="\t".join(map(str,row))+"\t"+ancestral+'\n'
    f_out.write(line)
    
f_out.close()

