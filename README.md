# Local PCAngsd

[![Documentation](https://img.shields.io/badge/-Documentation-blue)](https://alxsimon.github.io/local_pcangsd)
[![License](https://img.shields.io/github/license/alxsimon/local_pcangsd)](https://github.com/alxsimon/local_pcangsd/LICENCE.txt)

The objective of this python module is to combine [PCAngsd](https://github.com/Rosemeis/pcangsd) with [lostruct](https://github.com/jguhlin/lostruct-py)
to analyse low-coverage data (genotype likelihoods) with local PCA.

Genotype likelihood files can be large and will often not fit into memory.
This module leverages the use of [Xarray](https://xarray.dev/) to store and access genotype likelihoods on disk, in a data structure comparable to [`sgkit`](https://pystatgen.github.io/sgkit/latest/).

Similarly, PCA results are stored as an xarray dataset for easy manipulation and storage.

## Install

Requirements can be seen in the [`conda_env.yaml`](https://github.com/alxsimon/local_pcangsd/blob/main/conda_env.yaml) file.

The easiest way to install local_pcangsd and its dependencies is through conda:

```bash
mamba env create -f conda_env.yaml
# OR conda env create -f conda_env.yaml
conda activate local_pcangsd
git clone https://github.com/alxsimon/local_pcangsd.git
pip install ./local_pcangsd
```

## Usage

If you want to add the conda environment as a jupyter kernel

```bash
conda activate local_pcangsd
python -m ipykernel install --user --name local_pcangsd
```

Example code is presented in [`example.ipynb`](https://github.com/alxsimon/local_pcangsd/blob/main/example.ipynb).
