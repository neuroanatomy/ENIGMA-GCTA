#!/usr/bin/env python3

"""Estimate heritability for datasets with various parameters"""

import os
import sys
import shutil
import math
from inspect import getsourcefile
from configsimu import MYPLINK, GCTB, GRM, PHE_DIR, HSQ, N_SNP, N_IND, N_ITER, PCS, USE_SBATCH, HSQ_DIR, PRU_DIR, NBPROC

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))))
import preprocessing  # noqa: E402
sys.path.pop(0)


def main():
    """Entry point if called as an executable"""

    if USE_SBATCH:
        mode = "sbatch"
    else:
        mode = "direct"

    gen_file = os.path.join(PRU_DIR, 'all')

    # estimate hsq with GCTA
    for hsq in HSQ:
        for n_snp in N_SNP:
            for n_ind in N_IND:
                phe_dir = os.path.join(PHE_DIR, "hsq_" + str(hsq) + "-snp_" + str(n_snp) + "-ind_" +
                                       str(n_ind))

                phe_file = os.path.join(phe_dir, "phe_\\i")
                ind_file = os.path.join(phe_dir, "ind_\\i.txt")

                print("Remove header of individual lists...")
                for i in range(1, N_ITER+1):
                    ind_file_i = ind_file.replace("\\i", str(i))
                    phe_file_i = phe_file.replace("\\i", str(i))

                    source_file = open(ind_file_i, 'r')
                    source_file.readline()
                    target_file = open(phe_file_i + ".ind", 'w')

                    shutil.copyfileobj(source_file, target_file)
                    source_file.close()
                    target_file.close()

                out_dir_hsq = os.path.join(HSQ_DIR + "_gctb", os.path.basename(phe_dir))

                if os.path.exists(out_dir_hsq):
                    shutil.rmtree(out_dir_hsq)
                os.makedirs(out_dir_hsq)

                print("estimating hsq...")
                sys.stdout.flush()

                # estimate required memory
                nsnp = sum(1 for _ in open(gen_file + ".bim"))
                mem = 4 * nsnp * n_ind + nsnp * 500
                mem_str = str(math.ceil(mem / 1e6)) + "M"

                preprocessing.run([MYPLINK, 'mpirun', '-np', str(NBPROC), '--oversubscribe', GCTB,
                                   "--bfile", gen_file,
                                   "--pheno", phe_file + ".phen",
                                   "--keep", phe_file + ".ind",
                                   "--qcovar", os.path.join(PHE_DIR, "age.txt"),
                                   "--qcovar", PCS,
                                   "--covar", os.path.join(PHE_DIR, "centre.txt"),
                                   "--covar", os.path.join(PHE_DIR, "sex.txt"),
                                   "--bayes", "S",
                                   "--out", os.path.join(out_dir_hsq, os.path.basename(phe_file))],
                                  mode=mode,
                                  slurm_par=["-J", "simu_hsq",
                                             "--qos", "normal",
                                             "-p", "common",
                                             "--mem", mem_str,
                                             "-c", str(NBPROC),
                                             "-D", out_dir_hsq,
                                             "-W"],
                                  array=range(1, N_ITER+1),
                                  check=False)
                print("creating zip file...")
                sys.stdout.flush()
                shutil.make_archive(out_dir_hsq, "zip", os.path.dirname(out_dir_hsq), os.path.basename(out_dir_hsq))


if __name__ == '__main__':
    main()
