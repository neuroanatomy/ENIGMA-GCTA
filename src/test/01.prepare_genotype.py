#!/usr/bin/env python3

"""Prepare test genotype"""

import os
import subprocess
# from glob import glob
import pandas as pd
import numpy as np
from configtest import GEN_DIR, RAW_DIR, PLINK, N_SNP


def main():
    """Entry point if called as an executable"""
    np.random.seed(2017)

    # extract individuals with european ancestry and convert bcf file to plink format
    panel_file = os.path.join(RAW_DIR, "1000genomes",
                              "integrated_call_samples_v3.20130502.ALL.panel")
    panel_table = pd.read_table(panel_file)
    eur_ind = panel_table.loc[panel_table.super_pop == "EUR", "sample"]
    keep_file = os.path.join(GEN_DIR, "eur.ind")
    keep_table = pd.DataFrame({'FID': eur_ind, 'IID': eur_ind})
    if not os.path.exists(GEN_DIR):
        os.makedirs(GEN_DIR)
    keep_table.to_csv(keep_file, sep='\t', index=False)
    bcf_file = os.path.join(RAW_DIR, "1000genomes",
                            "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf")
    plink_file = os.path.join(GEN_DIR, "1000genomes_eur")
    if not os.path.exists(plink_file + ".bed"):
        subprocess.run([PLINK, "--bcf", bcf_file, "--keep", keep_file, "--make-bed", "--out",
                        plink_file, "--allow-extra-chr", "--maf", "0.005"], check=True)

    # update plink fam file from 1000 genomes ped file
    fam_table = pd.read_table(plink_file + ".fam", sep=' ',
                              names=['Family ID', 'Individual ID', 'Paternal ID', 'Maternal ID',
                                     'Gender', 'Phenotype'])
    fam_table.set_index('Individual ID', drop=False, inplace=True)
    ped_file = os.path.join(RAW_DIR, "1000genomes",
                            "integrated_call_samples_v2.20130502.ALL.ped")
    ped_table = pd.read_table(ped_file, index_col=1)
    fam_table.update(ped_table)
    if not os.path.exists(plink_file + ".fam.bak"):
        os.rename(plink_file + ".fam", plink_file + ".fam.bak")
    fam_table.to_csv(plink_file + ".fam", sep=' ', na_rep=-9, index=False, header=False)

    # extract SNP subset
    snp_table = pd.read_table(plink_file + ".bim", na_values=".", header=None)
    snp_table.dropna(inplace=True)
    snp_list = snp_table.sample(N_SNP).iloc[:, 1]
    snp_file = os.path.join(GEN_DIR, "mysnps.txt")
    snp_list.to_csv(snp_file, index=False)
    out_file = os.path.join(GEN_DIR, "all")
    subprocess.run([PLINK, "--bfile", plink_file, "--extract", snp_file, "--make-bed",
                    "--out", out_file], check=True)

    # remove intermediate plink file to free space
    # for file in glob(plink_file + ".*"):
    #     os.remove(file)


if __name__ == '__main__':
    main()
