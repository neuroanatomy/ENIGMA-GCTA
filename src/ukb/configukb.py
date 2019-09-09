"""Make config"""

import subprocess
import os

def getannexdir():
    """call git executable to find annex dir"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return subprocess.run(["git", "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE,
                          cwd=script_dir).stdout.decode().strip()


ANNEX_DIR = globals().get('ANNEX_DIR', getannexdir())
PLINK = globals().get('PLINK', 'plink')
PLINK2 = globals().get('PLINK2', 'plink2')
GCTA = globals().get('GCTA', 'gcta64')


DATASET_ID = '24220'

# indicate which dataset to use
DATASET = 'ukb-' + DATASET_ID

# data subdirectories
RAW_DIR = os.path.join(ANNEX_DIR, 'data', 'raw', )
DERIVED_DIR = os.path.join(ANNEX_DIR, 'data', 'derived')
EXTERNAL_DIR = os.path.join(ANNEX_DIR, 'data', 'external')

# directory containing all phenotypes
PHE_DIR = os.path.join(DERIVED_DIR, DATASET, '00.phenotype')

# derived data subdirectories
GEN_DIR = os.path.join(DERIVED_DIR, DATASET, '01.genotype')

# binaries directory
BIN_DIR = os.path.join(ANNEX_DIR, 'bin', 'ukb')

# UKB phenotype data
UKB_DATA = os.path.join(RAW_DIR, "ukb" + DATASET_ID + "-phenotype", "ukb" + DATASET_ID + ".enc")


# number of processors
NBPROC = 4

# remove list for individuals
# REMOVE_IND = os.path.join(RAW_DIR, "ukb" + DATASET_ID + "-phenotype", "w18584_20180503.csv")
REMOVE_IND = None

REGIONS = {
    "sex": "f.31.0.0",
    "height": "f.50.2.0",
    "centre": "f.54.2.0",
    "intelligence": "f.20016.2.0",
    "ICV": "f.25000.2.0",
    "brain": "f.25010.2.0",
    "thalamus": ["f.25011.2.0", "f.25012.2.0"],
    "caudate": ["f.25013.2.0", "f.25014.2.0"],
    "putamen": ["f.25015.2.0", "f.25016.2.0"],
    "pallidum": ["f.25017.2.0", "f.25018.2.0"],
    "hippocampus": ["f.25019.2.0", "f.25020.2.0"],
    "amygdala": ["f.25021.2.0", "f.25022.2.0"],
    "accumbens": ["f.25023.2.0", "f.25024.2.0"],
    "date.init": "f.53.0.0",
    "age.init": "f.21003.0.0",
    "date.mri": "f.53.2.0",
    "age.mri": "f.21003.2.0"
}

REGIONS_FREESURFER = {
    "ICV": ["aseg_eTIV"],
    "brain": ["aseg_BrainSegVol"],
    "thalamus": ["aseg_Left-Thalamus-Proper_Volume_mm3", "aseg_Right-Thalamus-Proper_Volume_mm3"],
    "caudate": ["aseg_Left-Caudate_Volume_mm3", "aseg_Right-Caudate_Volume_mm3"],
    "putamen": ["aseg_Left-Putamen_Volume_mm3", "aseg_Right-Putamen_Volume_mm3"],
    "pallidum": ["aseg_Left-Pallidum_Volume_mm3", "aseg_Right-Pallidum_Volume_mm3"],
    "hippocampus": ["aseg_Left-Hippocampus_Volume_mm3", "aseg_Right-Hippocampus_Volume_mm3"],
    "amygdala": ["aseg_Left-Amygdala_Volume_mm3", "aseg_Right-Amygdala_Volume_mm3"],
    "accumbens": ["aseg_Left-Accumbens-area_Volume_mm3", "aseg_Right-Accumbens-area_Volume_mm3"],
}

LOG10 = ["accumbens", "amygdala", "brain", "caudate", "height", "hippocampus", "ICV", "pallidum",
         "putamen", "thalamus"]
