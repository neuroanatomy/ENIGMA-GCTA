[![CircleCI](https://circleci.com/gh/neuroanatomy/genomic-architecture/tree/master.svg?style=svg)](https://circleci.com/gh/neuroanatomy/genomic-architecture/tree/master)

# Genomic architecture

This repository contains the code used for the analyses described in our manuscript "Polygenic architecture of human neuroanatomical diversity" (Preprint: https://www.biorxiv.org/content/10.1101/592337v3).

We used this code to analyse the genomic architecture of neuroanatomical diversity using magnetic resonance imaging and SNP data from > 26,000 individuals.


# Project organisation
The project organisation is inspired by that of the [Data Science Cookie Cutter](http://drivendata.github.io/cookiecutter-data-science). Raw data is stored in `/data/raw` and converted depending on its particularities to a dataset that is stored in `/data/derived`. Here, the dataset `/data/derived/test` is provided as illustration and for testing purposes.

# Instructions

## Load git submodules
```
git submodule init
git submodule update
```

## Install bioconda
```
conda config --add channels bioconda
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels defaults
```

## Create a conda environment (recommended)
```
conda create --name genomic-architecture python=3
conda activate genomic-architecture
```

## Install requirements
```
conda install --file requirements.txt
conda install --file requirements_r.txt
```

## Run analyses for test dataset
```
python3 src/00.run_all.py src/config/test.yml
```
