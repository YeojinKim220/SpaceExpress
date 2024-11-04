# SpaceExpress
SpaceExpress

## News
2024.05.07  

## Overview

## Getting started
See [Tutorials](./docs/source/notebook/).

For training one sample, check Tutorial 1. 

For training multiple sample, check Tutorial 3. 

## Software dependencies
Check requirement

## Installation


```bash
conda create -n spaceexpress python=3.11 scikit-learn pandas matplotlib jupyter conda-forge::scanpy conda-forge::python-igraph r::rpy2 conda-forge::r-lmtest conda-forge::r-fitdistrplus cconda-forge::r-dplyr -y

conda activate spaceexpress
pip3 install torch

cd SpaceExpress

python -m pip install -e .
```

<!-- python setup.py build

python setup.py install -->

## Citation
