# SpaceExpress
SpaceExpress

## News
2024.05.07  

## Overview

## Getting started
See [Tutorials](./docs/source/notebook/).
For training one sample, check [Tutorial 1](./docs/source/notebook/Tutorial_1.ipynb). 

## Software dependencies
Check requirement

## Installation
```bash
conda create -n spaceexpress python=3.11 scikit-learn pandas matplotlib jupyter conda-forge::scanpy conda-forge::python-igraph r::rpy2 conda-forge::r-lmtest conda-forge::r-fitdistrplus cconda-forge::r-dplyr -y
pip3 install torch

conda env create -f environment.yml
conda activate spaceexpress-env

cd SpaceExpress
python -m pip install -e .
```

<!-- python setup.py build

python setup.py install -->

## Citation
