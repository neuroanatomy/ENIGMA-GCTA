#!/usr/bin/env python3

"""Compute GWASs
"""

import os
import sys
import subprocess
import numpy as np


def getannexdir():
    """call git executable to find annex dir"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return subprocess.run(["git", "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE,
                          cwd=script_dir).stdout.decode().strip()


sys.path.insert(0, os.path.join(getannexdir(), 'src'))


import preprocessing
import config_dataset


def main(config_file):
    """Entry point if called as an executable"""

    config = config_dataset.config_dataset(config_file)

    plink2 = os.path.join(config.annex_dir, 'bin/plink2_linux_avx2_20190527/plink2')
    memory = "10000"

    in_prefix = 'all'
    in_file_allsnps = os.path.join(config.fil_dir, in_prefix)
    out_dir_gwas_allsnps = os.path.join(config.gwa_dir, 'gwas-all-sans')
    log_dir_gwas_allsnps = os.path.join(out_dir_gwas_allsnps, 'log')
    in_file_prunedsnps = os.path.join(config.pru_dir, in_prefix)
    out_dir_gwas_prunedsnps = os.path.join(config.gwa_dir, 'gwas-pruned-sans')
    log_dir_gwas_prunedsnps = os.path.join(out_dir_gwas_prunedsnps, 'log')
    pcs = os.path.join(config.phe_dir, "PC.txt")

    # use sex, center and age as covariates

    related_file = os.path.join(config.phe_dir, "ukb_all-0.025.rem")
    related_ind = np.loadtxt(related_file, usecols=[1], dtype=str)
    all_ind = np.loadtxt(in_file_allsnps + ".fam", usecols=[1], dtype=str)
    ukb_22110_list = os.path.join(config.derived_dir, "ukb-22110_maf001", "01.genotype", "all.fam")
    ukb_22110_ind = np.loadtxt(ukb_22110_list, usecols=[1], dtype=str)
    ukb_24220_list = os.path.join(config.derived_dir, "ukb-24220_maf001", "01.genotype", "all.fam")
    ukb_24220_ind = np.loadtxt(ukb_24220_list, usecols=[1], dtype=str)

    os.makedirs(log_dir_gwas_allsnps, exist_ok=True)
    os.makedirs(log_dir_gwas_prunedsnps, exist_ok=True)

    keep_ind = np.setdiff1d(all_ind, related_ind)
    remove_ind = np.setdiff1d(ukb_24220_ind, ukb_22110_ind)
    keep_ind = np.setdiff1d(keep_ind, remove_ind)
    keep_file = os.path.join(out_dir_gwas_allsnps, "keep.txt")
    np.savetxt(keep_file, np.column_stack((keep_ind, keep_ind)), fmt="%s")

    # slurm configuration
    if config.use_sbatch:
        smode = "sbatch"
    else:
        smode = "direct"

    for pheno in config.phe_list:

        # All SNPs
        out_prefix = 'all.' + pheno
        preprocessing.run([config.myplink, plink2,
                           "--bfile", in_file_allsnps,
                           "--keep", keep_file,
                           "--allow-no-sex", "--linear", "hide-covar",
                           "--pheno", os.path.join(config.phe_dir, pheno+".txt"),
                           "--qcovar", os.path.join(config.phe_dir, "age.txt"),
                           "--qcovar", pcs,
                           "--covar", os.path.join(config.phe_dir, "centre.txt"),
                           "--qcovar", os.path.join(config.phe_dir, "sex.txt"),
                           "--covar-variance-standardize",
                           "--ci", str(0.95),
                           "--threads", str(config.nbproc),
                           "--memory", memory,
                           "--out", os.path.join(out_dir_gwas_allsnps, out_prefix)],
                          mode=smode,
                          slurm_par=["-J", "gwas",
                                     "-D", log_dir_gwas_allsnps,
                                     "-c", str(config.nbproc),
                                     "--mem", memory + "M",
                                     "-p", "ghfc",
                                     "--qos", "ghfc"])

        # Pruned SNPs
        out_prefix = 'all.' + pheno
        preprocessing.run([config.myplink, plink2,
                           "--bfile", in_file_prunedsnps,
                           "--keep", keep_file,
                           "--allow-no-sex", "--linear", "hide-covar",
                           "--pheno", os.path.join(config.phe_dir, pheno+".txt"),
                           "--qcovar", os.path.join(config.phe_dir, "age.txt"),
                           "--qcovar", pcs,
                           "--covar", os.path.join(config.phe_dir, "centre.txt"),
                           "--qcovar", os.path.join(config.phe_dir, "sex.txt"),
                           "--covar-variance-standardize",
                           "--ci", str(0.95),
                           "--threads", str(config.nbproc),
                           "--memory", memory,
                           "--out", os.path.join(out_dir_gwas_prunedsnps, out_prefix)],
                          mode=smode,
                          slurm_par=["-J", "gwas",
                                     "-D", log_dir_gwas_prunedsnps,
                                     "-c", str(config.nbproc),
                                     "--mem", memory + "M",
                                     "-p", "ghfc",
                                     "--qos", "ghfc"])


if __name__ == '__main__':
    main(sys.argv[1])
