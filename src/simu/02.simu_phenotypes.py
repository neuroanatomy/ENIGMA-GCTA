#!/usr/bin/env python3

"""Simulate phenotypes for datasets with various parameters"""

import os
import sys
import shutil
from inspect import getsourcefile
import numpy as np
import pandas as pd
from configsimu import MYGCTA, GCTA, GRM, PRU_DIR, PHE_DIR, HSQ, N_SNP, N_IND, N_ITER, PCS, USE_SBATCH, HSQ_DIR

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))))
import preprocessing  # noqa: E402
sys.path.pop(0)


def main():
    """Entry point if called as an executable"""
    # initiate seed for pandas (but not for GCTA)
    np.random.seed(2017)

    gen_file = os.path.join(PRU_DIR, 'all')

    if USE_SBATCH:
        mode = "sbatch"
    else:
        mode = "direct"

    if not os.path.exists(PHE_DIR):
        os.makedirs(PHE_DIR)

    # load SNP table
    snp_table = pd.read_table(gen_file + ".bim", na_values=".", header=None).iloc[:, 1]
    snp_table.dropna(inplace=True)

    # load individual table
    ind_table = pd.read_table(PCS, sep=' ', dtype={'FID': str, 'IID': str}).loc[:, ['FID', 'IID']]

    # simulate heritable phenotypes with GCTA
    for hsq in HSQ:
        for n_snp in N_SNP:
            for n_ind in N_IND:
                out_dir = os.path.join(PHE_DIR, "hsq_" + str(hsq) + "-snp_" + str(n_snp) + "-ind_" +
                                       str(n_ind))
                if os.path.exists(out_dir):
                    shutil.rmtree(out_dir)
                os.makedirs(out_dir)

                phe_file = os.path.join(out_dir, "phe_\\i")
                snp_file = os.path.join(out_dir, "snp_\\i.txt")
                ind_file = os.path.join(out_dir, "ind_\\i.txt")

                print("Generating SNP and individual lists...")
                for i in range(1, N_ITER+1):
                    # extract SNP subset
                    snp_file_i = snp_file.replace("\\i", str(i))
                    snp_list = snp_table.sample(n_snp)
                    snp_list.to_csv(snp_file_i, index=False)

                    # extract individual subset
                    ind_file_i = ind_file.replace("\\i", str(i))
                    ind_list = ind_table.sample(n_ind).sort_values(by=["FID", "IID"])
                    ind_list.to_csv(ind_file_i, index=False, sep='\t')

                print("Simulating phenotypes...")
                sys.stdout.flush()
                # simulate phenotypes with GCTA
                preprocessing.run([GCTA,
                                   "--bfile", gen_file,
                                   "--keep", ind_file,
                                   "--simu-qt",
                                   "--simu-hsq", str(hsq),
                                   "--simu-causal-loci", snp_file,
                                   "--out", phe_file],
                                  mode=mode,
                                  slurm_par=["-J", "simu_pheno",
                                             "--qos", "fast",
                                             "-p", "common,dedicated",
                                             "--mem", "4G",
                                             "-D", out_dir,
                                             "-W"],
                                  array=range(1, N_ITER+1))

                # print("creating zip file...")
                # sys.stdout.flush()
                # shutil.make_archive(out_dir, "zip", os.path.dirname(out_dir),
                #                     os.path.basename(out_dir))


if __name__ == '__main__':
    main()
