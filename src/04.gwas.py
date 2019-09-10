#!/usr/bin/env python3

"""Compute GWASs
"""

import os
import sys
import pandas as pd
# from configall import PLINK, MYPLINK, USE_SBATCH
import preprocessing
import config_dataset


def main(config_file):
    """Entry point if called as an executable"""

    config = config_dataset.config_dataset(config_file)

    in_prefix = 'all'
    in_file_allsnps = os.path.join(config.fil_dir, in_prefix)
    out_dir_gwas_allsnps = os.path.join(config.gwa_dir, 'gwas-all')
    log_dir_gwas_allsnps = os.path.join(out_dir_gwas_allsnps, 'log')
    in_file_prunedsnps = os.path.join(config.pru_dir, in_prefix)
    out_dir_gwas_prunedsnps = os.path.join(config.gwa_dir, 'gwas-pruned')
    log_dir_gwas_prunedsnps = os.path.join(out_dir_gwas_prunedsnps, 'log')

    # use sex, center and age as covariates

    # Create dummy coded centre table
    # filter individuals in centre.txt keeping only those in the genotype file
    fam_table = pd.read_table(in_file_allsnps + ".fam", delim_whitespace=True, dtype=str,
                              names=['FID', 'IID', 'PID', 'MID', 'Gender', 'Phenotype'])
    centre_table = pd.read_table(os.path.join(config.phe_dir, "centre.txt"), delim_whitespace=True,
                                 index_col=False, dtype=str)
    centre_table = centre_table[centre_table.IID.isin(fam_table.IID)]
    centre_table.to_csv(os.path.join(config.phe_dir, "centre.cov"), sep='\t', index=False)

    # slurm configuration
    if config.use_sbatch:
        smode = "sbatch"
    else:
        smode = "direct"

    # All SNPs

    os.makedirs(log_dir_gwas_allsnps, exist_ok=True)

    for pheno in config.phe_list:
        out_prefix = 'all.' + pheno
        preprocessing.run([config.myplink, config.plink,
                           "--bfile", in_file_allsnps,
                           "--allow-no-sex", "--linear", "hide-covar",
                           "--pheno", os.path.join(config.phe_dir, pheno+".txt"),
                           "--qcovar", os.path.join(config.phe_dir, "age.txt"),
                           "--qcovar", config.pcs,
                           "--covar", os.path.join(config.phe_dir, "centre.cov"),
                           "--qcovar", os.path.join(config.phe_dir, "sex.txt"),
                           "--ci", str(0.95),
                           "--out", os.path.join(out_dir_gwas_allsnps, out_prefix)],
                          mode=smode,
                          slurm_par=["-J", "gwas",
                                     "-D", log_dir_gwas_allsnps])

    # Pruned SNPs

    os.makedirs(log_dir_gwas_prunedsnps, exist_ok=True)

    for pheno in config.phe_list:
        out_prefix = 'all.' + pheno
        preprocessing.run([config.myplink, config.plink,
                           "--bfile", in_file_prunedsnps,
                           "--allow-no-sex", "--linear", "hide-covar",
                           "--pheno", os.path.join(config.phe_dir, pheno+".txt"),
                           "--qcovar", os.path.join(config.phe_dir, "age.txt"),
                           "--qcovar", config.pcs,
                           "--covar", os.path.join(config.phe_dir, "centre.cov"),
                           "--qcovar", os.path.join(config.phe_dir, "sex.txt"),
                           "--ci", str(0.95),
                           "--out", os.path.join(out_dir_gwas_prunedsnps, out_prefix)],
                          mode=smode,
                          slurm_par=["-J", "gwas",
                                     "-D", log_dir_gwas_prunedsnps])


if __name__ == '__main__':
    main(sys.argv[1])
