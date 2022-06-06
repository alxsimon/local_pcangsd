# TO DO

* switch to GNU v3 (cf Nayuki code)

# Local PCAngsd

In development...

The final objective is to combine [PCAngsd](https://github.com/Rosemeis/pcangsd) with [lostruct](https://github.com/jguhlin/lostruct-py)
to analyse low-coverage data (genotype likelihoods) with local PCA.

## Install

Requirements can be seen in the [`conda_env.yaml`](https://github.com/alxsimon/local_pcangsd/blob/main/conda_env.yaml) file.

The PCAngsd package is installed from a forked version that packages it for easier access to internal functions.
See https://github.com/alxsimon/pcangsd.

The easiest way to install local_pcangsd and its dependencies is through conda:

```bash
conda env create -f conda_env.yaml
# OR
mamba env create -f conda_env.yaml
conda activate local_pcangsd
git clone https://github.com/alxsimon/local_pcangsd.git
pip install ./local_pcangsd
```

## Usage

If you want to add the conda environment as a jupyter kernel

```bash
conda activate local_pcangsd
mamba install -y ipykernel
python -m ipykernel install --user --name local_pcangsd
```

Example code is presented in [`example.ipynb`](https://github.com/alxsimon/local_pcangsd/blob/main/example.ipynb).


## Enclosing circle code

Code produced by [@Nayuki](https://github.com/nayuki),
available [here](https://github.com/nayuki/Nayuki-web-published-code/blob/master/smallest-enclosing-circle/smallestenclosingcircle.py).