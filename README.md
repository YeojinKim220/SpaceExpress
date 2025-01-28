# SpaceExpress

SpaceExpress is a tool for differential spatial expression (DSE) analysis using intrinsic coordinate systems of tissues. It leverages spline models and likelihood ratio testing to identify significant spatial gene expression changes across different biological conditions or experimental groups.

## News
**2024.11.26**  
- Initial release of SpaceExpress, enabling robust differential spatial expression analysis with support for single and multiple replicates.

## Getting Started
To learn how to use SpaceExpress for your spatial transcriptomics data, follow the tutorials:

- See **[Tutorials](./docs/source/notebook/)**.
- **Single Replicate Analysis**: Check [Tutorial 1](./docs/source/notebook/Tutorial_1.ipynb).
- **Multiple Replicate Analysis**: Check [Tutorial 2](./docs/source/notebook/Tutorial_2.ipynb).

These tutorials will guide you through data preparation, model setup, and interpreting the results of differential spatial expression analysis.

## Software Dependencies
To ensure compatibility and optimal performance, please check the required software dependencies specified in the `environment.yml` file. These include key libraries such as `scanpy`, `pandas`, `numpy`, and `rpy2` for seamless integration with R-based statistical models.

## Installation
Install the SpaceExpress environment and package with the following commands:

```bash
git clone https://github.com/YeojinKim220/SpaceExpress.git
conda env create -f environment.yml
conda activate spaceexpress-env

cd SpaceExpress
python -m pip install -e .
```

This will create and activate a Conda environment with all necessary dependencies and install SpaceExpress in editable mode for local development.

<details>
  <summary><b>Note:</b> If the installation via the provided <code>environment.yml</code> file fails, you can manually set up the environment using the steps below:</summary>

  ```bash
  # Step 1: Create and activate a new Conda environment
  conda create -n spaceexpress-env python=3.11 scikit-learn pandas matplotlib jupyter scanpy rpy2 -y
  conda activate spaceexpress-env

  # Step 2: Install additional dependencies
  conda install -y -c conda-forge python-igraph r-lmtest r-fitdistrplus r-dplyr r-lme4

  # Step 3: Install Python packages using pip
  pip install pygam
  pip install torch
  ```
</details>

## Citation
If you use SpaceExpress in your research, please cite:
```
@article{kim2024spaceexpress,
  title={SpaceExpress: a method for comparative spatial transcriptomics based on intrinsic coordinate systems of tissues},
  author={Kim, Yeojin and Ojha, Abhishek and Schrader, Alex and Lee, Juyeon and Wu, Zijun and Traniello, Ian M and Robinson, Gene E and Han, Hee Sun and Zhao, Sihai D and Sinha, Saurabh},
  journal={bioRxiv},
  pages={2024--12},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```

