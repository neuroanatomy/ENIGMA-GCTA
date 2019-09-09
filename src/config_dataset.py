"""Make config
General configuration of paths and binaries
"""

# from configall import ANNEX_DIR
import os
import sys
import subprocess
import re
import yaml
import preprocessing #pylint: disable=import-error


class Config(object):
    """Configuration for ukb pipeline"""
    def __init__(self, config_dict):

        # parameters imported from yaml file
        self.dataset = config_dict['dataset']
        self.annex_dir = config_dict.get('annex_dir') or getannexdir()
        self.bayesrv2 = config_dict.get('bayesRv2') or 'bayesRv2'
        self.use_sbatch = config_dict.get('use_sbatch') or False
        self.reml_call = config_dict.get('reml_call') or '--reml'
        self.reml_bivar_call = config_dict.get('reml_bivar_call') or '--reml-bivar'
        self.nbproc = config_dict.get('nbproc') or 1
        self.phe_list = config_dict['phe_list']

        # binaries
        bin_dir = os.path.join(self.annex_dir, "bin")
        bin_subdirs = next(os.walk(bin_dir))[1]
        if sys.platform == "linux":
            plink_re = re.compile(r"plink_linux_x86_64_[0-9]")
            gcta_re = re.compile(r"gcta(?!.*_[a-z]*$)")
            gcta_path = "gcta64"
            self.prsice_bin = os.path.join(bin_dir, "PRSice", "PRSice_linux")

        elif sys.platform == "darwin":
            plink_re = re.compile(r"plink_mac_[0-9]")
            gcta_re = re.compile(r"gcta.*_mac")
            gcta_path = os.path.join("bin", "gcta64")
            self.prsice_bin = os.path.join(bin_dir, "PRSice", "PRSice_mac")
        else:
            raise Exception("Unknown platform: " + sys.platform)

        plink_dirs = sorted(filter(plink_re.match, bin_subdirs))
        if plink_dirs:
            plink = os.path.join(bin_dir, plink_dirs[-1], "plink")
        else:
            plink = 'plink'
        gcta_dirs = sorted(filter(gcta_re.match, bin_subdirs))
        if gcta_dirs:
            gcta = os.path.join(bin_dir, gcta_dirs[-1], gcta_path)
        else:
            gcta = 'gcta64'
        self.plink = config_dict.get('plink') or plink
        self.gcta = config_dict.get('gcta') or gcta
        self.prsice_r = os.path.join(bin_dir, "PRSice", "PRSice.R")
        self.prsice = "Rscript " + self.prsice_r + " --prsice " + self.prsice_bin

        # maf intervals
        self.maf_bounds = config_dict.get('maf_bounds') or [0.001, 0.05, 0.20, 0.35]
        self.maf = self.maf_bounds[0]
        self.maf_intervals = list(zip(self.maf_bounds[:-1], self.maf_bounds[1:]))
        if self.maf_bounds[-1] != 0.5:
            self.maf_intervals.append((self.maf_bounds[-1], .5))

        # wrapper for plink
        self.myplink = os.path.join(self.annex_dir, "bin", "myplink", "myplink.sh")

        # wrapper for gcta
        self.mygcta = os.path.join(self.annex_dir, "bin", "mygcta", "mygcta.sh")

        # wrapper for bayesR
        # self.mybayesr = os.path.join(self.annex_dir, "bin", "mybayesR", "mybayesR.sh")

        # data subdirectories
        self.raw_dir = os.path.join(self.annex_dir, 'data', 'raw')
        self.derived_dir = os.path.join(self.annex_dir, 'data', 'derived')
        self.external_dir = os.path.join(self.annex_dir, 'data', 'external')

        # print python version
        print(sys.version)

        # derived data subdirectories
        self.dataset_dir = os.path.join(self.derived_dir, self.dataset)
        self.gen_dir = os.path.join(self.dataset_dir, '01.genotype')
        self.fil_dir = os.path.join(self.dataset_dir, '02.filtered')
        self.pru_dir = os.path.join(self.dataset_dir, '03.pruned')
        self.grm_dir = os.path.join(self.dataset_dir, '04.grm')
        self.gwa_dir = os.path.join(self.dataset_dir, '05.gwas')
        self.hsq_dir = os.path.join(self.dataset_dir, '06.hsq')

        # permutations
        if config_dict.get('permu_dir'):
            self.permu_dir = os.path.join(self.derived_dir, self.dataset, config_dict.get('permu_dir'))
        else:
            self.permu_dir = os.path.join(self.hsq_dir, 'permutations')
        self.permu_type = config_dict.get('permu_type') or 'simple'
        self.clean_permu = config_dict.get('clean_permu') or True
        self.nbit = config_dict.get('nbit') or 100
        self.array_lim = config_dict.get('array_lim') or 20

        # GRM cutoff
        self.grm_cutoff = 0.025

        # directory containing all phenotypes
        self.phe_dir = os.path.join(self.derived_dir, self.dataset, '00.phenotype')

        # subjects filtered out after QC filtering, used during variant pruning and filtering
        if os.path.exists(os.path.join(self.phe_dir, 'outliers.txt')):
            self.exclude_ind = os.path.join(self.phe_dir, 'outliers.txt')
        else:
            self.exclude_ind = None

        # SNPs range to be excluded in analyses (because they are in long range LD region)
        self.exclude_range = os.path.join(self.raw_dir, "exclude_range.txt")
        # self.exclude_range = None

        # subjects to keep for heritability : those selected as unrelated (used for the --keep parameter of gcta)
        self.keep_ind = os.path.join(self.grm_dir, 'grm-all-' + str(self.grm_cutoff), 'all-' + str(self.grm_cutoff)
                                     + '.grm.id')
        # Principal componants to be included as covariates in GWAS and heritability analyses
        self.pcs = os.path.join(self.grm_dir, 'grm-all-' + str(self.grm_cutoff), 'all-' + str(self.grm_cutoff)
                                + '.pca.eigenvec')

        # gene lists
        self.selgenes_neurodev = os.path.join(self.external_dir, 'neurodev.txt')
        self.selgenes_cnsexpression = os.path.join(self.external_dir, 'cnsexpression.txt')

        # ids of the phenotypes of interest, dir_pheno/id.txt will be used to find the phenotype file.
        # PHE_LIST = ["ICV", "brain", "accumbens", "amygdala", "caudate", "hippocampus",
        #            "pallidum", "putamen", "thalamus",  "height", "intelligence"]

        # self.phe_list = ["ICVlog10", "brainlog10", "accumbenslog10", "amygdalalog10", "caudatelog10",
        # "hippocampuslog10", "pallidumlog10", "putamenlog10", "thalamuslog10", "height", "intelligence"]

        # Phenotype and covariables
        # Check that files exist and contain enough non missing value

        # minimal number of non missing values accepted, arbitrary number for now
        self.min_nonmissing = 99

        self.phe_list = config_dict['phe_list']
        self.phe_list = [phenofile for phenofile in self.phe_list
                         if os.path.lexists(os.path.join(self.phe_dir, phenofile + '.txt'))
                         and preprocessing.pheno_count(os.path.join(self.phe_dir, phenofile + '.txt'))[0]
                         > self.min_nonmissing]

        self.quant_covar = [os.path.join(self.phe_dir, 'age.txt'), self.pcs]
        self.quant_covar = [covarfile for covarfile in self.quant_covar
                            if os.path.lexists(covarfile)
                            and preprocessing.pheno_count(covarfile)[0] > self.min_nonmissing]

        self.qual_covar = [os.path.join(self.phe_dir, 'sex.txt'), os.path.join(self.phe_dir, 'centre.txt')]
        self.qual_covar = [covarfile for covarfile in self.qual_covar
                           if os.path.lexists(covarfile)
                           and preprocessing.pheno_count(covarfile)[0] > self.min_nonmissing]


def getannexdir():
    """call git executable to find annex dir"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return subprocess.run(["git", "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE,
                          cwd=script_dir).stdout.decode().strip()


def config_dataset(config_file):
    """Generate config from yaml file.
    Required fields are:
      - dataset

    :param str config_file: yaml configuration file
    :return: config object
    :rtype: Config
    """

    config_dict = yaml.safe_load(open(config_file))

    config = Config(config_dict)

    print("Phenotypes used for dataset " + config.dataset + " are: " + ', '.join(config.phe_list))
    print("Qualitative covariables used for dataset " + config.dataset + " are: " + ', '.join(config.qual_covar))
    print("Quantitative covariables used for dataset " + config.dataset + " are: " + ', '.join(config.quant_covar))

    return config
