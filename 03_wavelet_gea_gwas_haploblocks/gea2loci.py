#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEA outlier SNPs of Todesco et al. (2020), run on mrfox.
For each taxon and each climate or soil variable, SNPs with |Bayes factor| > 20 in the BayPass runs without
inversion covariates are written as one BED per variable (gea_selection/<variable>_<taxon>.bed).
These files are copied to the local gea_loci/gea_selection/ folder and read by gea_genes.R.
"""

sps=["annuus","argophyllus","petfal","petpet"]
import pandas as pd
import os
import csv
def tell_sig_gea(sp,bf,header):
	bf.columns=header.columns
	gene_list=list()
	for col in range(3,header.shape[1]):
		out_file=f"/TEMP_SCRATCH/zhiqin/selective_sweeps/enrich/gea_selection/{bf.columns[col]}_{sp}.bed"
		if not os.path.exists(out_file):	
			aim=bf[abs(bf.iloc[:,col])>20]
			if not aim.empty:   
				out=aim.iloc[:,0].str.split("_",expand=True)
				out['start']=out.iloc[:,2].astype(int)-1
				out.iloc[:,[0,3,2]].to_csv(out_file,sep='\t',index=False,header=False)

for sp in sps: 
	df_file_dir="/moonriseNFS/wild_gwas/gea_results/climate/snps_HA412/without_inv_covariates/"+sp
	bf_file=df_file_dir+"/var_out_"+sp+'_2018_HA412_maf_0.03_all_covariates_BayPass_IS_STD_noinv.txt.gz'
	bf=pd.read_csv(bf_file,sep="\t",header=None)
	header=pd.read_csv(f"/moonriseNFS/wild_gwas/gea_results/climate/snps_HA412/with_inv_covariates/{sp}/column.names.{sp}.varout.txt",sep="\t")
	tell_sig_gea(sp,bf,header)
###for soil
	so_file_dir="/moonriseNFS/wild_gwas/gea_results/soil/snps_HA412/without_inv_covariates/"+sp
	bf_file=f"{so_file_dir}/{sp}_HA412_invfree_matrix_all_covariates_baypass.BF.txt.gz"
	bf=pd.read_csv(bf_file,sep=" ",header=None)
	header=pd.read_csv("/moonriseNFS/wild_gwas/gea_results/soil/header_info.txt",sep="\t")
	tell_sig_gea(sp,bf,header)
