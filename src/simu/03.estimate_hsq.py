#!/usr/bin/env python3

"""Estimate heritability for datasets with various parameters"""

import os
import sys
import shutil
from inspect import getsourcefile
from configsimu import MYGCTA, GCTA, GRM, PHE_DIR, HSQ, N_SNP, N_IND, N_ITER, PCS, USE_SBATCH, HSQ_DIR

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))))
import preprocessing  # noqa: E402
sys.path.pop(0)


def main():
    """Entry point if called as an executable"""

    if USE_SBATCH:
        mode = "sbatch"
    else:
        mode = "direct"

    if not os.path.exists(PHE_DIR):
        os.makedirs(PHE_DIR)

    # estimate hsq with GCTA
    for hsq in HSQ:
        for n_snp in N_SNP:
            for n_ind in N_IND:
                phe_dir = os.path.join(PHE_DIR, "hsq_" + str(hsq) + "-snp_" + str(n_snp) + "-ind_" +
                                       str(n_ind))

                phe_file = os.path.join(phe_dir, "phe_\\i")

                out_dir_hsq = os.path.join(HSQ_DIR, os.path.basename(phe_dir))

                if os.path.exists(out_dir_hsq):
                    shutil.rmtree(out_dir_hsq)
                os.makedirs(out_dir_hsq)

                print("estimating hsq...")
                sys.stdout.flush()
                preprocessing.run([MYGCTA, GCTA,
                                   "--grm-bin", GRM,
                                   "--pheno", phe_file + ".phen",
                                   "--qcovar", os.path.join(PHE_DIR, "age.txt"),
                                   "--qcovar", PCS,
                                   "--covar", os.path.join(PHE_DIR, "centre.txt"),
                                   "--covar", os.path.join(PHE_DIR, "sex.txt"),
                                   "--out", os.path.join(out_dir_hsq, os.path.basename(phe_file)),
                                   "--reml-no-constrain"],
                                  mode=mode,
                                  slurm_par=["-J", "simu_hsq",
                                             "--qos", "fast",
                                             "-p", "common,dedicated",
                                             "--mem", "2G",
                                             "-D", out_dir_hsq,
                                             "-W"],
                                  array=range(1, N_ITER+1),
                                  check=False)
                print("creating zip file...")
                sys.stdout.flush()
                shutil.make_archive(out_dir_hsq, "zip", os.path.dirname(out_dir_hsq), os.path.basename(out_dir_hsq))


if __name__ == '__main__':
    main()
