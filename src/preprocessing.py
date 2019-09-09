#!/usr/bin/env python3

"""PLINK preprocessing.

 .. _Google Python Style Guide:
   http://google.github.io/styleguide/pyguide.html
"""

import subprocess
import os
import sys
import math
import time
from struct import unpack, calcsize
import collections
import pandas as pd

# http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html


def main(in_dir, in_prefix,
         plink='plink'):
    """Entry point if called as an executable

    Args:
        in_dir (str): The folder containing the input data (.bed + .bim + .fam files).
        in_prefix (str): The prefix for the input data
        plink (str):      path of the PLINK software to be used
        in_dir/in_prefix.bed + in_dir/in_prefix.bim + in_dir/in_prefix.fam files
        are required by the PLINK commands.
    """

    filter_prune(in_dir, in_prefix, plink)


def linecount(filename):
    """count the number of lines of a file"""
    return int(os.popen('wc -l '+filename).read().split()[0])


def read_grm_bin_n(prefix):
    """read the number of SNPs used to compute a GRM from GCTA-format binary file
    adapted from an R function on GCTA website
    Davis McCarthy, June 2014.
    Copied from https://github.com/davismcc/GCTAtools/blob/master/GCTAtools.py
    """
    n_file_name = prefix + ".grm.N.bin"
    entry_format = 'f'  # N is stored as a float in the binary file
    entry_size = calcsize(entry_format)
    # Open the binary file
    with open(n_file_name, mode='rb') as file:
        # entry_count = os.fstat(f.fileno()).st_size / entry_size
        record = file.read(entry_size)
        num = unpack(entry_format, record)[0]
        num = int(num)
    return num


def pheno_count(pheno_file):

    """Read phenotype file and count occurences

    Args:
        pheno_file (str): path of phenotype file to be read

    Returns:
        list: first element is number of *non* missing values,
              second element contains the dictionary output of
              collections.Counter function that contains the number of
              occurences of each level of phenotype (-9 is for missing value).

    """

    dtf = pd.read_csv(pheno_file, sep='\\s+|\\t+', engine='python')
    annot = dtf.iloc[:, 2].values
    counter = collections.Counter(annot)
    # number of missing values: -9
    nb_missing = counter[-9]
    nb_nonmissing = len(annot) - nb_missing

    return [nb_nonmissing, counter]


def filter_prune(in_dir, in_prefix,
                 plink='plink',
                 maf=0.001,
                 other_plink_par=None):
    """Filter and prune SNPs

    Args:
        in_dir (str): The folder containing the input data (.bed + .bim + .fam files).
        in_prefix (str): The prefix for the input data
        plink (str):      path of the PLINK software to be used
        maf (float):    The minimum minor allele frequecy to be keeped
        other_plink_par (list): list of plink parameters

        in_dir/in_prefix.bed + in_dir/in_prefix.bim + in_dir/in_prefix.fam files
        are required by the PLINK commands.
    """

    # out_dir_filtering, out_prefix_filtering, out_dir_pruning, out_prefix_pruning):
    # https://github.com/chapmanb/bcbb/blob/master/nextgen/scripts/plink_to_vcf.py
    # base_dir = os.getcwd()

    # subprocess.call(["module load plink/1.90b2m"])

    # exp = os.path.basename(in_prefix)
    # exp = os.path.splitext(base)[0]

    gen_dir = os.path.join(in_dir, '01.genotype')
    fil_dir = os.path.join(in_dir, '02.filtered')
    pru_dir = os.path.join(in_dir, '03.pruned')

    gen_dirprefix = os.path.join(gen_dir, in_prefix)

    for ext in ['.fam', '.bed', '.bim']:
        if not os.path.exists(gen_dirprefix + ext):
            print('Input file ' + gen_dirprefix + ext + ' is missing')
            sys.exit(1)

    print('Output directory for filtering: ' + fil_dir)  # RT out_dir_filtering)
    print('Output directory for pruning: ' + pru_dir, flush=True)  # RT out_dir_pruning)

    if not os.path.exists(fil_dir):  # RT out_dir_filtering)
        os.makedirs(fil_dir)  # RT out_dir_filtering)

    if not os.path.exists(pru_dir):  # RT out_dir_pruning)
        os.makedirs(pru_dir)  # RT out_dir_pruning)

    # variant selection
    plink_variant_selection(bfile=os.path.join(gen_dir, in_prefix),  # RT in_dirprefix,
                            out_dir=fil_dir,  # RT out_dir_filtering,
                            out_prefix='all',
                            maf=maf, geno=0.05, hwe=1e-6, mind=0.1,
                            plink=plink,
                            other_plink_par=other_plink_par)
    # variant pruning
    plink_variant_pruning(bfile=os.path.join(fil_dir, in_prefix),
                          out_dir=pru_dir,  # RT out_dir_pruning,
                          out_prefix='all',
                          maf=maf,
                          geno=0.01,
                          hwe=1e-6,
                          mind=0.1,
                          indep_window=50,
                          indep_count=10,
                          indep_variance=10,
                          plink=plink,
                          other_plink_par=other_plink_par)

    # RT return out_prefix_pruning


def plink_variant_selection(bfile, out_dir, out_prefix='',
                            maf=0.001,
                            geno=0.05,
                            hwe=1e-6,
                            mind=0.1,
                            plink='plink',
                            other_plink_par=None):
    """Run PLINK for variant selection.

    Variant selection using as default the parameters:
        MAF > 1%,
        genotyping rate > 95%,
        significance threshold for HW equilibrium test > 1e-6,
        subject missing genotypes < 10%.

    Writes filtered .bed + .bim + .fam.

    Args:
        bfile (str):      The --bfile flag causes the binary fileset bfile.bed + bfile.bim
                          + bfile.fam to be referenced.
        out_dir (str):    where to writes the output files out.bed and out.bim  and
                          out.fam.
        out_prefix (str): prefix for the output files
        maf (float):      minimum allele frequency (default 0.01)
        geno (float):     filters out all variants with missing call rates exceeding the
                          provided value (default 0.05) to be removed.
        hwe (float):      filters out all variants which have Hardy-Weinberg equilibrium
                          exact test p-value below the provided threshold (default 1e-6).
                          We recommend setting a low threshold-serious genotyping errors
                          often yield extreme p-values like 1e-50 which are detected by
                          any reasonable configuration of this test, while genuine
                          SNP-trait associations can be expected to deviate slightly from
                          Hardy-Weinberg equilibrium (so it's dangerous to choose a
                          threshold that filters out too many variants).
        mind (float):     filters out all samples with missing call rates exceeding the
                          provided value (0.1).
        plink (str):      path of the PLINK software to be used
        other_plink_par (list): list of plink parameters

        Other args used in call to PLINK:
            --noweb
            --make-bed creates a new PLINK 1 binary fileset,
                after applying sample/variant filters and other operations below.
                For example. {--bfile}.bed + .bim + .fam.


    Returns:
        boolean:        True if successful, raise error otherwise.

    Raises:
        OSError:        Error when running plink.
    """

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out_base = os.path.join(out_dir, out_prefix)

    if other_plink_par is None:
        other_plink_par = []

    if not os.path.exists("{0}.bed".format(out_base)):  # Do not overwrite
        command = [plink,
                   "--bfile", str(bfile),
                   "--noweb",
                   "--maf", str(maf),
                   "--geno", str(geno),
                   "--hwe", str(hwe),
                   "--mind", str(mind),
                   "--autosome",
                   "--make-bed",
                   "--out", str(out_base)
                   ] + other_plink_par
        print(" ".join(command), flush=True)
        try:
            subprocess.check_call(command)
        except OSError as ex:
            raise ex
        return True


def plink_make_set(bfile, out_file,
                   gene_set,
                   subset,
                   border=0,
                   plink='plink'):
    """Run PLINK for selection of variants.

    Variant selection using as default the parameters:
        MAF > 1%,
        genotyping rate > 95%,
        significance threshold for HW equilibrium test > 1e-6,
        subject missing genotypes < 10%.

    Writes filtered .bed + .bim + .fam.

    Args:
        bfile (str): The --bfile flag causes the binary fileset bfile.bed + bfile.bim + bfile.fam
          to be referenced.
        out_file (str): prefix for the output files
        gene_set (str): file listing the genomic (genic) regions to be considered
        border (int): border/margin to considered outside the genomic regions
        subset (str): file containing subset of gene ids
        plink (str): path of the PLINK software to be used

        Other args used in call to PLINK:
            --write-snplist

    Returns:
        boolean: True if successful, raise error otherwise.

    Raises:
        OSError:        Error when running plink.
    """

    out_dir = os.path.dirname(out_file)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    command = [plink,
               "--bfile", str(bfile),
               "--noweb",
               "--make-set", str(gene_set),
               "--make-set-border", str(border),
               "--subset", str(subset),
               "--gene-all",            # remove also variants not in gene_set (not a default in plink 1.9)
               "--write-snplist",
               "--out", str(out_file)
               ]
    print(" ".join(command), flush=True)
    try:
        subprocess.check_call(command)
    except OSError as ex:
        raise ex
    return True


def plink_variant_pruning(bfile, out_dir, out_prefix='',
                          maf=0.001,
                          geno=0.01,
                          hwe=1e-6,
                          mind=0.1,
                          indep_window=50,
                          indep_count=5,
                          indep_variance=10,
                          plink='plink',
                          other_plink_par=None):

    """Run PLINK for variant filtering/pruning.

    Variant selection using as default the parameters:
        MAF > 1%,
        genotyping rate > 95%,
        significance threshold for HW equilibrium test > 1e-6,
        subject missing genotypes < 10%
    Exclusion of variants in strong linkage desiquilibirum (R^2 > 0.9) within a window of 50 SNPs.

    It produces a pruned subset of markers that are in approximate linkage equilibrium
    with each other, writing the IDs to *.prune.in (and the IDs of all excluded variants
    to *.prune.out). They are currently based on correlations between genotype allele counts;
    phase is not considered.

    Writes filtered .bed + .bim + .fam.

    To do:
        add argument keep later
        --keep ../../../10imagen2089/data/to-include.txt
        accepts a space/tab-delimited text file with family IDs in the first column
        and within-family IDs in the second column,
        and removes all unlisted samples from the current analysis.

    Args:
        bfile (str): The --bfile flag causes the binary fileset bfile.bed + bfile.bim + bfile.fam
          to be referenced.
        out_dir (str): directory for output files
        out_prefix (str): prefix for output files
        maf (float): minimum allele frequency (default 0.05)
        geno (float): filters out all variants with missing call rates exceeding the provided value
        (default 0.01) to be removed.
        hwe (float): filters out all variants which have Hardy-Weinberg equilibrium exact test
          p-value below the provided threshold (default 1e-6). We recommend setting a low threshold
          as serious genotyping errors often yield extreme p-values like 1e-50 which are detected
          by any reasonable configuration of this test, while genuine SNP-trait associations can
          be expected to deviate slightly from Hardy-Weinberg equilibrium (so it's dangerous to
          choose a threshold that filters out too many variants).
        mind (float): filters out all samples with missing call rates exceeding the provided
          value (0.1).
        indep PLINK argument that requires three parameters:
        indep_window (Union[int, str]): a window size in variant count units (default 50),
        indep_count (Union[int, str]): a variant count to shift the window at the end of each step
            (default 5),
        indep_variance (int): a variance inflation factor (VIF) threshold (default 10).
            At each step, all variants in the current window with VIF exceeding the threshold
              re removed.
        plink (str): path of the PLINK software to be used
        other_plink_par (list): list of plink parameters

        Other args used in call to PLINK:
            --noweb
            --make-bed creates a new PLINK 1 binary fileset,
                after applying sample/variant filters and other operations below.
                For example. {--bfile}.bed + .bim + .fam.

    Returns:
        boolean: True if successful, raise error otherwise.

    Raises:
        OSError: when plink raise an error
    """

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out_base = os.path.join(out_dir, out_prefix)

    if other_plink_par is None:
        other_plink_par = []

    if not os.path.exists("{0}.bed".format(out_base)):  # Do not overwrite

        niter = 3

        for i in range(1, niter+1):

            if i > 1:
                extract_par = ["--extract", out_base+".pruning{}.prune.in".format(i-1)]
            else:
                extract_par = []

            try:
                # variant pruning, this writes .in and .out files
                subprocess.check_call([plink,
                                       "--bfile", bfile,
                                       "--noweb",
                                       "--maf", str(maf),
                                       "--geno", str(geno),
                                       "--hwe", str(hwe),
                                       "--mind", str(mind),
                                       "--indep", str(indep_window),
                                       str(indep_count), str(indep_variance),
                                       "--out", out_base + ".pruning" + str(i)] + extract_par + other_plink_par)
            except OSError as ex:
                raise ex

        try:
            # use .in file to make bed files containing the selected variants
            subprocess.check_call([plink,
                                   "--bfile", bfile,
                                   "--noweb",
                                   "--make-bed",
                                   "--extract", out_base+'.pruning{}.prune.in'.format(niter),
                                   "--out", out_base] + other_plink_par)
        except OSError as ex:
            raise ex

        try:
            # Compute frequencies for the "pruned" data
            subprocess.check_call([plink,
                                   "--bfile", out_base,
                                   "--freq",
                                   "--noweb",
                                   "--allow-no-sex",  # By default, unless the input is loaded
                                   # with --no-sex1, samples with ambiguous sex have their
                                   # phenotypes set to missing when analysis commands are run.
                                   # Use --allow-no-sex to prevent this.
                                   "--out", out_base + ".freq"] + other_plink_par)
        except OSError as ex:
            raise ex
            # add argument keep later
            # --keep ../../../10imagen2089/data/to-include.txt
            # accepts a space/tab-delimited text file with family IDs in the first column
            # and within-family IDs in the second column,
            # and removes all unlisted samples from the current analysis.
    return True


def gcta_grm_pca(grm,
                 out,
                 nbpcs=10,
                 gcta='gcta'):
    """Run GCTA for PCA computation.

    Use GCTA to compute Principal Component Analysis of a GRM matrix:

    From the three output files (grm.bin, grm.N.bin, grm.id) generated by the function gcta_grm,
    output the first n (n = 20, by default) eigenvalues (saved as *.eigenval, plain text file)
    and eigenvectors (saved as *.eigenvec, plain text file).

    See http://cnsgenomics.com/software/gcta/pca.html for more details.

    Writes two output files (.eigenvec, .eigenval).

    Args:
        grm (str): path to the three 'grm' files (grm.bin, grm.N.bin, grm.id) generated by GCTA
        out (str): prefix for output files
        nbpcs (int): number of principcal components to estimate
        gcta (str): path of the GCTA software to be used



    Returns:
        boolean: True if successful, Raise error otherwise.

    Raises:
        OSError: when GCTA raises an error
    """

    try:
        subprocess.check_call([gcta,
                               "--grm-bin", str(grm),
                               "--pca", str(nbpcs),
                               "--out", str(out)])
    except OSError as ex:
        raise ex
    return True

# def gcta_exec(in_file,
#               out_file,
#               par_input='--bfile',
#               gcta='gcta',
#               ncpus=1,
#               other_gcta_par=None,
#               sbatch=False):

#     """Run GCTA"""

#     if other_gcta_par is None:
#         other_gcta_par = []

#     out_dir = os.path.dirname(out_file)

#     if not os.path.exists(out_dir):
#         os.makedirs(out_dir)

#     log_dir = os.path.join(out_dir, "log")

#     if sbatch:
#         bmode = "sbatch"
#         if not os.path.exists(log_dir):
#             os.makedirs(log_dir)
#     else:
#         bmode = "direct"

#     try:
#         args = [gcta,
#                 par_input, in_file,
#                 "--out", out_file,
#                 "--thread-num", str(ncpus)] + other_gcta_par
#         slurm_par = ["-J", "gcta_grm",
#                      "-p", "dedicated,common",
#                      "-D", log_dir,
#                      "-c", str(ncpus)]
#         jid = run(args, mode=bmode, slurm_par=slurm_par)
#         return jid
#     except OSError as ex:
#         raise ex


def gcta_grm_part(in_file,
                  out_file,
                  par_input='--bfile',
                  gcta='gcta',
                  ncpus=1,
                  other_gcta_par=None,
                  sbatch=False,
                  nparts=None):
    """Compute GRM in multiple parts"""

    def file_len(fname):
        """Compute number of lines of a file"""
        num_lines = sum(1 for _ in open(fname))
        return num_lines

    if other_gcta_par is None:
        other_gcta_par = []

    out_dir = os.path.dirname(out_file)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    log_dir = os.path.join(out_dir, "log")

    if sbatch:
        bmode1 = "sbatch"
        bmode2 = "sbatch"
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
    else:
        bmode1 = "direct"
        bmode2 = "shell"

    max_mem = 5

    if par_input != "--bfile":
        raise ValueError("Error: unknown value '{}' for parameter 'par_input'".format(par_input))
    if nparts is None:
        # compute number of part for not exceeding 5 GB of memory per task
        nind = file_len(in_file + ".fam")
        nparts = math.ceil(nind**2 * 6 / (1024**3 * max_mem * 0.9))

    if nparts == 1:
        print("Make GRM once", flush=True)
        try:
            args = [gcta,
                    par_input, in_file,
                    "--make-grm",
                    "--out", out_file,
                    "--thread-num", str(ncpus)
                    ] + other_gcta_par
            slurm_par = ["-J", "gcta_grm",
                         "-p", "common",
                         "--qos", "normal",
                         "-D", log_dir,
                         "-c", str(ncpus)]
            jid = run(args, mode=bmode1, slurm_par=slurm_par)
            return jid
        except OSError as ex:
            raise ex
    else:
        print("Dividing GRM {} in {} parts".format(out_file, nparts), flush=True)
        try:
            args = [gcta,
                    par_input, in_file,
                    "--make-grm-part", str(nparts), "\\i",
                    "--out", out_file,
                    "--thread-num", str(ncpus)
                    ] + other_gcta_par
            slurm_par = ["-J", "gcta_grm",
                         "-p", "common",
                         "--qos", "normal",
                         "--mem", str(max_mem) + "G",
                         "-D", log_dir,
                         "-c", str(ncpus)]
            jid = run(args, mode=bmode1, slurm_par=slurm_par, array=range(1, nparts+1))
        except OSError as ex:
            raise ex
        cmd = ("cat "+out_file+".part_"+str(nparts)+"_*.grm.id > "+out_file+".grm.id;" +
               "cat "+out_file+".part_"+str(nparts)+"_*.grm.bin > "+out_file+".grm.bin;" +
               "cat "+out_file+".part_"+str(nparts)+"_*.grm.N.bin > "+out_file+".grm.N.bin;" +
               "rm "+out_file+".part_"+str(nparts)+"_*.grm.{id,bin,N.bin}")
        slurm_par = ["-J", "gcta_grm_merge",
                     "-p", "common",
                     "--qos", "normal",
                     "-D", log_dir,
                     "--dependency", "afterok:" + jid]
        jid = run(cmd, mode=bmode2, slurm_par=slurm_par)
    return jid


def gcta_grm(in_file,
             out_file,
             par_input='--bfile',
             gcta='gcta',
             per_chr=False,
             chr_list=None,
             ncpus=1,
             other_gcta_par=None,
             sbatch=False):

    """Use GCTA to compute genetic relationship matrix (GRM):

    From three bed/bim/fam test files, estimates the genetic relationship matrix (GRM)
    between pairs of individuals from a set of SNPs and save the lower triangle elements of
    the GRM to binary files.

    Writes three output files
    test.grm.bin (it is a binary file which contains the lower triangle elements of the GRM).
    test.grm.N.bin (it is a binary file which contains the number of SNPs used to calculate the
      GRM).
    test.grm.id (no header line; columns are family ID and individual ID, see above).
    see http://cnsgenomics.com/software/gcta/estimate_grm.html for more details.

    Args:
        in_file (str): depending on par_input, either bfile fileset referencing bfile.bed
            + bfile.bim + bfile.fam, or grm file set referencing .grm.bin .grm.id .grm.N.bin files
            from gcta.
        out_file (str): prefix for output files
        gcta (str): path of the GCTA software to be used
        per_chr (bool): True if the GRM have to be computed by chromosome first before they are
            merged, False is default. This option is very useful to deal with large dataset.
        chr_list (list): list of chromosomes for which a GRM should be computed
        other_gcta_par (list): list of additional gcta parameters to be used
        par_input (str): "--bfile" (default) if the GRM is computed from bfiles (bed/bim/fam),
            "--grm-bin" if it is computed from a previous GRM output
        ncpus (int): number of cpus to use
        sbatch (boolean): True to use sbatch to compute the GRM.

        Other args used in call to GCTA:
            --automosome to not use X or Y chromosomes in the estimation

    Returns:
        boolean: True if successful, raise error otherwise.

    Raises:
        OSError: When GCTA raises an error.
    """

    if chr_list is None:
        chr_list = list(range(1, 23))
    if other_gcta_par is None:
        other_gcta_par = []

    out_dir = os.path.dirname(out_file)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    log_dir = os.path.join(out_dir, "log")

    if sbatch:
        rmode = "srun"
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
    else:
        rmode = "direct"

    print('Output directory for GRM: ' + out_dir, flush=True)

    # Compute GRM
    if per_chr:
        # GRM by chromosome, then merge
        # file that will contain path of chromosome GRM files
        myfile = open(os.path.join(out_dir, 'all.multi.txt'), 'w')
        jids = []

        for chrom in chr_list:
            grmchr_file = out_file + '-chr' + str(chrom)
            if not os.path.isfile(grmchr_file + '.grm.bin'):
                try:
                    jids += [gcta_grm_part(in_file=in_file,
                                           out_file=grmchr_file,
                                           par_input=par_input,
                                           gcta=gcta,
                                           ncpus=ncpus,
                                           other_gcta_par=["--chr", str(chrom)] + other_gcta_par,
                                           sbatch=sbatch)]
                except OSError as ex:
                    raise ex
            myfile.write("%s\n" % grmchr_file)
        myfile.close()

        # merge GRM by chromosome
        if not os.path.isfile(out_file + '.grm.bin'):
            try:
                run([gcta,
                     "--mgrm-bin", os.path.join(out_dir, 'all.multi.txt'),
                     "--make-grm-bin",
                     "--out", out_file,
                     "--thread-num", str(ncpus)
                     ] + other_gcta_par,
                    mode=rmode,
                    slurm_par=["-J", "gcta_grm_merge",
                               "-p", "common,dedicated",
                               "--qos", "fast",
                               "-D", log_dir,
                               "-c", str(ncpus),
                               "--dependency", "afterok:" + ":".join(jids)])
            except OSError as ex:
                raise ex
    else:
        # GRM all chromosomes at the same time
        if not os.path.isfile(out_file + 'grm.bin'):
            jid = gcta_grm_part(in_file=in_file,
                                out_file=out_file,
                                par_input=par_input,
                                gcta=gcta,
                                ncpus=ncpus,
                                other_gcta_par=["--autosome"] + other_gcta_par,
                                sbatch=sbatch)
        run(["printf", "GRM computed.\n"],
            mode=rmode,
            slurm_par=["-J", "gcta_completed",
                       "-p", "common",
                       "--qos", "normal",
                       "--dependency", "afterok:" + jid])

    return True


def gcta_grm_filter(in_file,
                    out_file,
                    grm_cutoff=0.025,
                    gcta='gcta',
                    per_chr=False,
                    chr_list=None,
                    ncpus=1,
                    other_gcta_par=None,
                    srun=False):

    """Use GCTA to compute genetic relationship matrix (GRM):

    From a GRM file, filter individuals to keep only a subset with no pair of individuals having
        a relatedness over the value of grm_cutoff.

    Writes three output files
    test.grm.bin (it is a binary file which contains the lower triangle elements of the GRM).
    test.grm.N.bin (it is a binary file which contains the number of SNPs used to calculate the
      GRM).
    test.grm.id (no header line; columns are family ID and individual ID, see above).
    see http://cnsgenomics.com/software/gcta/estimate_grm.html for more details.

    Args:
        in_file (str): a grm file set referencing .grm.bin .grm.id .grm.N.bin files from gcta.
        out_file (str): prefix for output files
        grm_cutoff (float): if more than zero,  used for parameter --grm-cutoff to remove one of a
            pair of individuals with estimated relatedness larger than the specified cut-off value
            (e.g. 0.025). This is ran *after* the GRM for all SNPs is computed.
        gcta (str): path of the GCTA software to be used
        per_chr (bool): True if the GRM have been computed by chromosome first before they were
            merged, False is default.
        chr_list (list): list of chromosomes for which the GRM have been computed
        ncpus (int): number of cpus to use
        srun (boolean): if True, launch GCTA with srun.
        other_gcta_par (list): list of additional gcta parameters to be used

        Other args used in call to GCTA:
            --automosome to not use X or Y chromosomes in the estimation

    Returns:
        boolean: True if successful, raise error otherwise.

    Raises:
        OSError: raised if GCTA crashes
    """

    if chr_list is None:
        chr_list = list(range(1, 23))
    if other_gcta_par is None:
        other_gcta_par = []

    if srun:
        rmode = "srun"
    else:
        rmode = "direct"

    out_dir = os.path.dirname(out_file)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # filter the merged result for grm_cutoff

    print('Output directory for filtered GRM: ' + out_dir)
    sys.stdout.flush()

    try:
        run([gcta,
             "--grm-bin", in_file,
             "--grm-cutoff", str(grm_cutoff),
             "--make-grm-bin",
             "--out", out_file
             ] + other_gcta_par,
            mode=rmode,
            slurm_par=[ "-p", "common,dedicated",
                       "--qos", "fast"])
    except OSError as ex:
        raise ex

    if per_chr:
        # GRM by chromosome using filtered individuals
        keep_ind = out_file + '.grm.id'

        for chrom in chr_list:
            grmchr_file_filtered = out_file + '-chr' + str(chrom)
            grmchr_file = in_file + '-chr' + str(chrom)
            try:
                run([gcta,
                     "--grm-bin", grmchr_file,
                     "--keep", keep_ind,
                     "--make-grm-bin",
                     "--out", grmchr_file_filtered
                     ] + other_gcta_par,
                    mode=rmode,
                    slurm_par=[ "-p", "common,dedicated",
                               "--qos", "fast"])
            except OSError as ex:
                raise ex

    return True


def plink_exec(bfile, out_dir, out_prefix='',
               plink='plink',
               other_plink_par=None):
    """Run PLINK

        Args:
            bfile (str): The --bfile flag causes the binary fileset bfile.bed + bfile.bim
                + bfile.fam to be referenced.
            out_dir (str): where to writes the output files out.bed and out.bim  and out.fam.
            out_prefix (str): prefix for the output files
            plink (str): path of the PLINK software to be used
            other_plink_par (list): list of plink parameters

            Other args used in call to PLINK:
                --noweb

        Returns:
            boolean: True if successful, raise error otherwise.

        Raises:
            OSError: raised if PLINK fails to run
    """

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if other_plink_par is None:
        other_plink_par = []

    out_base = os.path.join(out_dir, out_prefix)
    print('Output file for plink: ' + out_base)
    sys.stdout.flush()

    plinkcall = [plink,
                 '--bfile', bfile,
                 '--out', out_base] + other_plink_par
    try:
        subprocess.check_call(plinkcall)
    except OSError as ex:
        raise ex
    return True


def gcta_hsq(in_file,
             out_file,
             par_input='--grm-bin',
             gcta='gcta',
             mygcta='mygcta',
             other_gcta_par=None,
             ncpus=1,
             sbatch=True,
             sbatch_par_j='hsq',
             sbatch_par_p='common,dedicated',
             sbatch_par_qos='fast',
             sbatch_par_mem=None,
             sbatch_par_gres='disk:100'):

    """Use GCTA to compute the heritability estimates:

    From three grm.bin/grm.N.bin/grm.id  files, estimates the proportion of variance
    in a phenotype explained by all SNPs (i.e. the SNP-based heritability);

    Writes a .hsq file containing the
    see http://cnsgenomics.com/software/gcta/#GREMLanalysis for more details.

    Args:
        in_file (str): grm file set referencing .grm.bin .grm.id .grm.N.bin files
            from gcta.
        out_file (str): prefix for output files
        par_input (str): "--grm-bin" if it is computed from one grm,
                         "--mgrm-bin" if it is computed from multiple grms
        gcta (str): path of the GCTA software to be used
        mygcta (str): path of mygcta wrapper
        other_gcta_par (list): list of additional gcta parameters to be used
        ncpus (int): number of cpus to use.
        sbatch (bool): True if gcta has to be run through sbatch, False if it is run in local
        sbatch_par_j (string): sbatch job name
        sbatch_par_p (str): sbatch partition
        sbatch_par_qos (string): sbatch qos
        sbatch_par_mem (string): sbatch allocated memory
        sbatch_par_gres (string): sbatch parameter for disk allocation


    Returns:
        bool: True if successful, raise error otherwise.

    Raises:
        OSError: raised if GCTA crashes
        ValueError: raised if par_input not reconized

    """

    if other_gcta_par is None:
        other_gcta_par = []

    out_dir = os.path.dirname(out_file)
    out_base = os.path.basename(out_file)
    log_dir = os.path.join(out_dir, "log")

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    print('Output directory for h^2: ' + out_dir)
    sys.stdout.flush()

    if sbatch_par_mem is None:
        # Estimate required memory in kilobytes
        # see docs/hsq-mem.txt for experimental measures on required memory
        if par_input == "--grm-bin":
            grm = in_file
            ncomp = 2
        elif par_input == '--mgrm-bin':
            with open(in_file) as fobj:
                grm = fobj.readline().strip()
                # ncomp = number of genetic components + one for the environmental component
                ncomp = 2 + sum(1 for _ in fobj)    # 2 because the first line is already read
        else:
            raise ValueError("Error: unknown value '{}' for parameter 'par_input'"
                             .format(par_input))
        nind = linecount(grm + ".grm.id")
        if '--reml-bivar' in other_gcta_par:
            ncomp *= 3
            # nind *= 2.2   # overestimate nind because reml-bivar seems to need more memory
            nind *= 2.2

        mem = math.ceil(500000 + (0.025 + 0.009 * ncomp) * nind**2)
        sbatch_par_mem = str(mem) + 'K'
        # do not use less than 1 cpu per 10GB
        ncpus = max(int(ncpus), int(math.ceil(mem/1e7)))

    if sbatch:
        smode = "sbatch"
    else:
        smode = "direct"

    try:
        run([mygcta, gcta,
             par_input, in_file,
             "--out", out_file,
             "--thread-num", str(ncpus)] + other_gcta_par,
            mode=smode,
            check=False,
            slurm_par=["-J", sbatch_par_j,
                       "-p", sbatch_par_p,
                       "--qos", sbatch_par_qos,
                       "--mem", sbatch_par_mem,
                       "--gres", sbatch_par_gres,
                       "-D", log_dir,
                       "-o", out_base + "-%j.out",
                       "-e", out_base + "-%j.out",
                       "-c", str(ncpus)])

    except OSError as ex:
        raise ex

    return True


# def submit_sbatch_gcta(in_file,
#                        out_file,
#                        gcta='gcta',
#                        mygcta='mygcta',
#                        sbatch_par_J='hsq',
#                        sbatch_par_p='dedicated -p common',
#                        sbatch_par_qos='fast',
#                        sbatch_par_mem='1G',
#                        sbatch_par_gres='disk:100',
#                        sbatch_par_cpuspertask = '1',
#                        gcta_par=''):

#     """Run GCTA through a sbatch call:

#     Using a preformatted string containing the input parameter for the GCTA software,
#     run the GCTA command on the cluster through sbatch.

#     Writes a result .hsq file and gcta log files in the output folder.
#     Writes a run_sbatch.out file in the current folder containing the sbatch stdout.

#     Args:
#         in_file (str): grm file set referencing .grm.bin .grm.id .grm.N.bin files
#             from gcta.
#         out_file (str): prefix for output files
#         gcta (str): path of the GCTA software to be used
#         gcta_par (str): string containing the call to gcta (only parameters)
#         sbatch_par_* (str): dictionary providing sbatch parameters -J --qos --mem --gres


#     Returns:
#         True if successful, raise error otherwise.
#     """


#     bn = os.path.basename(out_file)
#     dn = os.path.dirname(out_file)

# #     sbatch_command  = """sbatch  -J  hsq \
# # -p dedicated -p common \
# # --qos   fast \
# # --mem   1G \
# # --gres  disk:100 \
# # -e  {bn}.out \
# # -o  {bn}.out \
# # -D {dn}/log \
# # {callbin} \
# # {bin} """.format(bn = bn, dn = dn, callbin=MYGCTA, bin = gcta)

#     sbatch_command = """sbatch  -J  {J} \
# -p {p} \
# --qos   {qos} \
# --mem   {mem} \
# --gres  {gres} \
# --cpus-per-task {cpus} \
# -e  {bn}.out \
# -o  {bn}.out \
# -D {dn}/log \
# {callbin} \
# {bin} \
# {command} """.format(J = sbatch_par_J, p = sbatch_par_p, qos = sbatch_par_qos,
# mem =  sbatch_par_mem, gres = sbatch_par_gres, cpus = sbatch_par_cpuspertask, bn = bn,
# dn = dn, callbin=mygcta, bin = gcta, command =  gcta_par)


#     print(sbatch_command)
#     sys.stdout.flush()

#     # juste transcrire les boucles en python
#     # tester si la fonction submit marche

#     with open('run_sbatch.out', "a") as f:

#         runsbatch = subprocess.Popen(sbatch_command,
#                                      stdout=f,
#                                      stderr=f,
#                                      shell=True)

#         runsbatch.wait()
#         # capturer sbatch return code ($?) pour savoir si commande a ete bien lancee
#         if runsbatch.returncode != 0:
#             print('sbatch command', sbatch_command, 'failed.', file=sys.stderr)

def run(args,
        mode="direct",
        slurm_par=None,
        array=None,
        check=True):

    """Run program with SLURM or directly
    """

    if slurm_par is None:
        slurm_par = []

    if not isinstance(args, list):
        args = [args]

    if array is None:
        return run_single(args, mode, slurm_par, check=check)

    else:
        return run_array(args, mode, slurm_par, array, check=check)


def run_single(args, mode, slurm_par, check=True):
    """run single command"""
    if mode == "direct":
        subprocess.run(args, check=check)

    elif mode == "shell":
        subprocess.run(args, shell=True, check=check)

    elif mode == "srun":
        subprocess.run(["salloc"] + slurm_par + ["srun"] + args, check=check)

    elif mode == "sbatch":
        ret = subprocess.run(["sbatch"] + slurm_par + ["--wrap", " ".join(args)],
                             stdout=subprocess.PIPE, check=check)
        print(ret.stdout.decode().strip(), flush=True)
        job_id = ret.stdout.decode().rsplit()[3]
        return job_id

    else:
        raise ValueError("invalid value for mode parameter")

    return ""


def run_array(args, mode, slurm_par, array, array_limit=None, check=True):
    """run several commands in array"""
    if mode == "sbatch":
        cmd = " ".join(args).replace("\\i", '${SLURM_ARRAY_TASK_ID}')
        cmd = cmd.replace("\\j", '${SLURM_ARRAY_JOB_ID}')
        print(cmd, flush=True)
        if array_limit is None:
            array_limit_command = ''
        else:
            array_limit_command = '%' + str(array_limit)

        ret = subprocess.run(["sbatch", "-a", ",".join(map(str, array)) + array_limit_command]
                             + slurm_par + ["--wrap", cmd], stdout=subprocess.PIPE, check=check)

        print(ret.stdout.decode().strip(), flush=True)
        job_id = ret.stdout.decode().rsplit()[3]

    else:
        job_id = str(int(time.time()))
        for i in array:
            args_i = [arg.replace("\\i", str(i)).replace("\\j", job_id) for arg in args]

            if mode == "direct":
                subprocess.run(args_i, check=check)

            elif mode == "shell":
                subprocess.run(args_i, shell=True, check=check)

            elif mode == "srun":
                subprocess.run(["salloc"] + slurm_par + ["srun"] + args_i, check=check)

            else:
                raise ValueError("invalid value for mode parameter")

    return job_id
