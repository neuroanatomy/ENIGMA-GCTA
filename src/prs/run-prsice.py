#!/usr/bin/env python3

"""Estimate genome wide polygenic scores with PRSice from PLINK GWAS outputs.
"""

import sys
import subprocess
import os
import pandas as pd


def getannexdir():
    """call git executable to find annex dir"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return subprocess.run(["git", "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE,
                          cwd=script_dir).stdout.decode().strip()


sys.path.insert(0, os.path.join(getannexdir(), 'src'))

# pylint: disable=E0401,C0413
import preprocessing   # noqa: E0401
import config_dataset  # noqa: E0401


def main(config_file):
    """Entry point if called as an executable"""

    config = config_dataset.config_dataset(config_file)

    gwas_suffix = "all"
    gwas_dataset = "ukb-22110_imputed"
    gwas_dir = os.path.join(config.derived_dir, gwas_dataset, "05.gwas", "gwas-" + gwas_suffix)
    target_file = os.path.join(config.gen_dir + "", "all")
    result_dir = os.path.join(config.dataset_dir, "PRSice", gwas_suffix + "." + gwas_dataset)
    remove_list = os.path.join(config.derived_dir, "ukb-22110_maf001", "01.genotype", "all.fam")

    os.makedirs(result_dir, exist_ok=True)

    if config.use_sbatch:
        mode = "sbatch"
    else:
        mode = "shell"

    quant_covar_table = pd.concat([pd.read_csv(f, index_col=["FID", "IID"], sep=r'\s+', dtype=str)
                                   for f in config.quant_covar], axis=1).dropna()
    qual_covar_table = pd.concat([pd.read_csv(f, index_col=["FID", "IID"], sep=r'\s+', dtype=str)
                                  for f in config.qual_covar], axis=1).dropna()
    covar_table = pd.concat([quant_covar_table, qual_covar_table], axis=1).dropna()

    cov_file = os.path.join(result_dir, "covariates.cov")
    covar_table.to_csv(cov_file, sep='\t', index=True)

    for pheno in config.phe_list:

        gwas_file = os.path.join(gwas_dir, "all." + pheno + ".assoc.linear")
        if not os.path.exists(gwas_file):
            print("Skipping GWAS not found for pheno: " + pheno)
            continue
        out_dir = os.path.join(result_dir, pheno)
        log_dir = os.path.join(out_dir, "log")
        os.makedirs(log_dir, exist_ok=True)

        command = [config.prsice,
                   "--base", gwas_file,
                   "--target", target_file,
                   "--pheno-file", os.path.join(config.phe_dir, pheno + ".txt"),
                   "--cov-file", cov_file,
                   "--cov-col", ','.join(covar_table.columns),
                   "--cov-factor", ','.join(qual_covar_table.columns),
                   "--remove", remove_list,
                   "--thread", str(config.nbproc),
                   "--out", os.path.join(out_dir, pheno)]
        print(" ".join(command))
        preprocessing.run(" ".join(command),
                          mode=mode,
                          check=False,
                          slurm_par=["-J", "PRSice",
                                     "-p", "ghfc",
                                     "--qos", "ghfc",
                                     "--mem", "120G",
                                     "-D", log_dir,
                                     "-o", os.path.join(log_dir, pheno + "-%j.out"),
                                     "-e", os.path.join(log_dir, pheno + "-%j.out"),
                                     "-c", str(config.nbproc)])


if __name__ == '__main__':
    main(sys.argv[1])
