# Local PCAngsd

In development...

The final objective is to combine [PCAngsd](https://github.com/Rosemeis/pcangsd) with [lostruct](https://github.com/jguhlin/lostruct-py)
to analyse low-coverage data (genotype likelihoods) with local PCA.

For now the lostruct part is not done.

## Requirements

- PCAngsd (my fork for now, requires a manual install): https://github.com/alxsimon/pcangsd
- h5py
- numpy
- pandas
- (lostruct at some point)

## Install

1) `pip install git+https://github.com/alxsimon/pcangsd.git`
2) `pip install h5py numpy pandas`
3) `pip install git+https://github.com/alxsimon/local_pcangsd.git`

## Usage

```python
import local_pcangsd as lp
```

**Some examples go here...**
