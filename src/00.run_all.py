#!/usr/bin/env python3

""" Master script that runs in the appropriate order all the scripts necessary to
    reproduce the derived data.
    @Author Roberto Toro, 1 Octobre 2017
"""

# @todo Should exec using
# slurm --mem 10GB [file.py] |tee file.log.txt
# @todo test or remove from scripts the installation of R packages
# R packages can already be installed with `conda install --file requirements_r.txt`

# pylint: disable=C0103
# or else we have to put everything in a main function


import os
import sys
import subprocess
import rpy2.robjects as robjects
import config_dataset


# === configurations of paths === #
config_file = sys.argv[1]
config = config_dataset.config_dataset(config_file)

src_dir = os.path.join(config.annex_dir, "src")

subprocess.run([os.path.join(src_dir, "01-01.install_binaries.py"), config_file], check=True)
subprocess.run([os.path.join(src_dir, "01-02.download_external.py"), config_file], check=True)

# === Filtering and variant pruning === #
subprocess.run([os.path.join(src_dir, "02.filter_prune.py"), config_file], check=True)

# === Compute GRM for all variants and variant partitions === #
subprocess.run([os.path.join(src_dir, "03-01.grm_all.py"), config_file], check=True)
subprocess.run([os.path.join(src_dir, "03-02.grm_pca.py"), config_file], check=True)
subprocess.run([os.path.join(src_dir, "03-03.grm_genic.py"), config_file], check=True)
subprocess.run([os.path.join(src_dir, "03-04.grm_maf.py"), config_file], check=True)
subprocess.run([os.path.join(src_dir, "03-05.grm_gene_lists.py"), config_file], check=True)


# === GWAS === #
subprocess.run([os.path.join(src_dir, "04.gwas.py"), config_file], check=True)


# === Compute heritability estimates === #
subprocess.run([os.path.join(src_dir, "05-01.hsq_all.py"), config_file], check=True)
subprocess.run([os.path.join(src_dir, "05-02.hsq_nopca.py"), config_file], check=True)
subprocess.run([os.path.join(src_dir, "05-03.hsq_genic.py"), config_file], check=True)
subprocess.run([os.path.join(src_dir, "05-04.hsq_maf.py"), config_file], check=True)
subprocess.run([os.path.join(src_dir, "05-05.hsq_gene_lists.py"), config_file], check=True)
subprocess.run([os.path.join(src_dir, "05-06.hsq_bivariate.py"), config_file], check=True)
subprocess.run([os.path.join(src_dir, "05-07.hsq_perchr.py"), config_file], check=True)


# === Run permutations to estimate heritability enrichment in the different partitions === #
subprocess.run([os.path.join(src_dir, "05-08.hsq_permutations.py"), config_file], check=True)


# === Extract the hsq results and merge them in text files === #
hsqsummarydir = os.path.join(config.hsq_dir, 'hsq-summary')
os.makedirs(hsqsummarydir, exist_ok=True)
subprocess.run([os.path.join(src_dir, "06-01.summarize_hsq.sh"), config.hsq_dir, config.grm_dir,
                os.path.join(hsqsummarydir, 'summary')], check=True)


# === Plot GWAS results === #

r_source = robjects.r['source']

r_source(os.path.join(src_dir, "06-02.plot_gwas.R"))
plotGWAS = robjects.globalenv['plotGWAS']
plotGWAS(os.path.join(config.gwa_dir, 'gwas-all'))  # gwas-pruned ?

# === Plot heritability estimates === #
# plots in R using python package rpy2
# the plot scripts rely on the following R packages:
# data.table, ggplot2, foreach, GenomicRanges, ggbio

# to be tested
# Note: the required packages are added to requirements_r.txt, so maybe no need for these lines
# if install_Rpackages:
#     utils = rpackages.importr('utils')
#     base = rpackages.importr('base')
#     utils.chooseCRANmirror(ind=31) # Paris mirror
#     utils.chooseBioCmirror(ind=31)
#     # install CRAN packages
#     utils.install_packages("data.table")
#     utils.install_packages("ggplot2")
#     # install Bioconductor packages
#     base.source("http://www.bioconductor.org/biocLite.R")
#     bioclite = robjects.globalenv['biocLite']
#     biocLite('GenomicRanges')
#     biocLite('ggbio')

# ** for heritability results **
r_source = robjects.r['source']

r_source(os.path.join(src_dir, "06-03.plot_hsq.R"))
plotHSQ = robjects.globalenv['plotHSQ']
plotHSQ(prefix_hsq_summary_file=os.path.join(hsqsummarydir, 'summary_'),
        hsq_outdir=os.path.join(config.hsq_dir, ''))

# ** plots and writes univariate heritability results **
r_source(os.path.join(src_dir, "process_hsq_nopartition.R"))
readHsq = robjects.globalenv['readHsq']
plot_hsq = robjects.globalenv['plot_hsq']

#pandas2ri.activate()

resallsnps = readHsq(path_res=os.path.join(config.hsq_dir, 'hsq-all'))
resallsnps.to_csvfile(path=os.path.join(hsqsummarydir, 'hsq-all.txt'), sep='\t', quote=False, row_names=False)
plot_hsq(hsq=os.path.join(hsqsummarydir, 'hsq-all.txt'), pattern_pheno_ignore='left|right',
         outpdf=os.path.join(hsqsummarydir, 'hsq-all.pdf'))

resallsnps_nopca = readHsq(path_res=os.path.join(config.hsq_dir, 'hsq-nopca'))
resallsnps_nopca.to_csvfile(path=os.path.join(hsqsummarydir, 'hsq-nopca.txt'),
                            sep='\t', quote=False, row_names=False)
ggallsnps_nopca = plot_hsq(hsq=os.path.join(hsqsummarydir, 'hsq-nopca.txt'),
                           pattern_pheno_ignore='left|right',
                           outpdf=os.path.join(hsqsummarydir, 'hsq-nopca.pdf'))

# ** plots and writes heritability results partitioned by chromosome **
r_source(os.path.join(src_dir, "process_hsq_perchr.R"))
write_plot_hsq_perchr = robjects.globalenv['write_plot_hsq_perchr']
write_plot_hsq_perchr(hsq_outdir=config.hsq_dir, output_dir=hsqsummarydir)

# ** plots and writes partitioned heritability results **
r_source(os.path.join(src_dir, "process_partition.R"))
write_plot_hsq_partition = robjects.globalenv['write_plot_hsq_partition']
for type_partition, group_partition in zip(['maf', 'genic', 'cnsexpression', 'neurodev'],
                                           ['maf', 'updown-margin20-50', 'cnsexpression', 'neurodev']):

    res = write_plot_hsq_partition(hsq_outdir=config.hsq_dir, output_dir=hsqsummarydir,
                                   type_partition=type_partition, group_partition=group_partition)

# ** plots and writes genetic/phenotypic/environmental correlations**
r_source(os.path.join(src_dir, "process_bivariate.R"))
process_bivariate = robjects.globalenv['process_bivariate']
res = process_bivariate(path_res=os.path.join(config.hsq_dir, 'hsq-biv'),
                        prefix_res="all\\.rg\\\\?=0\\..*.hsq",
                        output_file=os.path.join(hsqsummarydir, 'summary_bivariate.txt'))
plot_bivariate = robjects.globalenv['plot_bivariate']
plot_bivariate(res=os.path.join(hsqsummarydir, 'summary_bivariate.txt'),
               dirout=hsqsummarydir)
