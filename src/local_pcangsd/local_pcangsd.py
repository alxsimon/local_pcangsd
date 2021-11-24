import numpy as np
import pandas as pd
import xarray as xr
import sgkit as sg
from sgkit.window import window_statistic
from pcangsd_core import shared
from pcangsd_core import covariance
import os, contextlib

DIM_VARIANT = "variants"
DIM_SAMPLE = "samples"
DIM_ALLELE = "alleles"
DIM_GENOTYPE = "genotypes"
DIM_PC = "PCs"
DIM_WINDOW = "windows"

def create_dataset(df: pd.DataFrame, contigs: list) -> xr.Dataset:
    """
    df: pandas DataFrame (chunk of the beagle file)
    contigs: list of contig names
    """
    sample_id = np.array(df.columns[range(3, df.shape[1], 3)])

    variants = df.marker.str.rsplit('_', 1, expand=True)
    variants.columns = ['contig', 'position']

    flat_gl = df.iloc[:,3:].to_numpy(dtype=np.float64)
    ncol = flat_gl.shape[1]
    GL = np.stack(
        [
            flat_gl[:, range(0, ncol, 3)],
            flat_gl[:, range(1, ncol, 3)],
            flat_gl[:, range(2, ncol, 3)],
        ],
        axis=-1
    )

    variant_contig = np.array([contigs.index(n) for n in variants.contig], dtype=np.int64)

    ds = xr.Dataset(
        data_vars={
            'variant_contig_name': ([DIM_VARIANT], variants.contig.to_numpy()),
            'variant_contig': ([DIM_VARIANT], variant_contig),
            'variant_position': ([DIM_VARIANT], variants.position.to_numpy(dtype=np.int64)),
            'sample_id': ([DIM_SAMPLE], sample_id),
            'allele': ([DIM_VARIANT, DIM_ALLELE], df.iloc[:,[1,2]].to_numpy(dtype=np.int8)),
            'genotype_likelihood': ([DIM_VARIANT, DIM_SAMPLE, DIM_GENOTYPE], GL),
        },
        attrs=dict(contigs=contigs),
    )
    return ds


def beagle_to_zarr(input: str, store: str, chunksize: int=10000) -> None:
    """
    Converts an ANGSD genotype likelihood dataset to a Zarr array on disk.
    """

    # first pass to obtain contig names
    df_chunked = pd.read_csv(input, sep='\t', chunksize=chunksize)
    contigs = [df.marker.str.rsplit('_', 1, expand=True).iloc[:,0].unique() for df in df_chunked]
    contigs = list(np.unique(np.hstack(contigs)))

    df_chunked = pd.read_csv(input, sep='\t', chunksize=chunksize)
    for i, df in enumerate(df_chunked):
        ds = create_dataset(df, contigs=contigs)
        if i == 0:
            ds.to_zarr(store, mode="w")
        else:
            ds.to_zarr(store, append_dim='variants')


def load_dataset(store: str, chunksize: int=10000) -> xr.Dataset:
    ds = xr.open_zarr(store, chunks=chunksize)
    return ds


def window(ds: xr.Dataset, type: str, size: int) -> xr.Dataset:
    """
    Wrapper arround sgkit.window_by_position.
    Size either in bp or number of variants depending on type.
    """

    if type not in {"position", "variant"}:
        raise ValueError(
            "Window type accepted are only 'position' or 'variant'."
        )
    
    if type == "position":
        return sg.window_by_position(ds, size=size)
    elif type == "variant":
        return sg.window_by_variant(ds, size=size)


def pcangsd_wrapper(
    gl: np.array, k: int,
    emMAF_iter: int=200, emMAF_tole: float=1e-4, emMAF_t: int=1, 
    emPCA_e: int=0, emPCA_iter: int=100, emPCA_tole: float=1e-5, emPCA_t: int=1,
) -> tuple[np.array, float, np.array, np.array]:
    """
    Wrapper around PCAngsd function returning a tuple
    similar to what is returned by lostruct.eigen_windows.

    gl: genotype likelihoods as formatted in the beagle_to_zarr output
        i.e. (variants, samples, genotypes=3)
    """

    L = gl[:,:,:-1].reshape(gl.shape[0], -1)
    L = L.astype('float32')
    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
            f = shared.emMAF(L, iter=emMAF_iter, tole=emMAF_tole, t=emMAF_t)
            C, P, _ = covariance.emPCA(L, f, e=emPCA_e, iter=emPCA_iter, tole=emPCA_tole, t=emPCA_t)
    vals, vectors = np.linalg.eig(C)
    total_variance = np.sum(np.power(C, 2).flatten())
    return C, total_variance, vals[:k], vectors[:, :k].T


def pca_windows(ds: xr.Dataset) -> np.array:
    """
    Run PCAngsd on each window.
    Return results as np.array of objects (lostruct style).
    """
    
    if 'window_start' not in ds or 'window_stop' not in ds:
        raise Exception("Variables 'window_start' and 'window_stop' not defined in the Dataset.")

    window_starts = ds.window_start.values
    window_stops = ds.window_stop.values
    result = []
    for start, stop in zip(window_starts, window_stops):
        result.append(pcangsd_wrapper(ds.genotype_likelihood[start:stop].values, k=10))
    
    return np.array(result, dtype=object) # same result as vstack


#===========================================================
# CODE NOT WORKING BELOW

def pcangsd_wrapper_vectors(gl: np.array) -> np.array:
    L = gl[:,:,:-1].reshape(gl.shape[0], -1)
    L = L.astype('float32')
    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
            f = shared.emMAF(L, iter=10, tole=1e-4, t=1)
            C, P, K = covariance.emPCA(L, f, e=0, iter=10, tole=1e-4, t=1)
    vals, vectors = np.linalg.eig(C)
    return vectors


def pca_windows_sgkit_stat(ds: xr.Dataset) -> xr.Dataset:
    """
    Run PCAngsd on each window.
    Only return the vectors (for now).
    """
    
    if 'window_start' not in ds or 'window_stop' not in ds:
        raise Exception("Variables 'window_start' and 'window_stop' not defined.")

    res = window_statistic(
        ds.genotype_likelihood,
        pcangsd_wrapper_vectors,
        window_starts=ds.window_start.values,
        window_stops=ds.window_stop.values,
        dtype=np.float32,
    )
    return res.compute(scheduler='single-threaded') # parallel version not working


#====================================================

# statistic function to apply to each window
def pca_gufunc(gl: np.array,
    emMAF_iter: int=200, emMAF_tole: float=1e-4, emMAF_t: int=1, 
    emPCA_e: int=0, emPCA_iter: int=100, emPCA_tole: float=1e-5, emPCA_t: int=1,
) -> tuple:
    L = gl[:,:,:-1].reshape(gl.shape[0], -1)
    L = L.astype('float32')

    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
            f = shared.emMAF(L, iter=emMAF_iter, tole=emMAF_tole, t=emMAF_t)
            C, P, K = covariance.emPCA(L, f, e=emPCA_e, iter=emPCA_iter, tole=emPCA_tole, t=emPCA_t)
    vals, vectors = np.linalg.eig(C)

    return (vectors, vals, C, P, K)

def pca_block(
    ds: xr.Dataset, window_size: int,
    emMAF_iter: int=200, emMAF_tole: float=1e-4, emMAF_t: int=1, 
    emPCA_e: int=0, emPCA_iter: int=100, emPCA_tole: float=1e-5, emPCA_t: int=1,
) -> xr.Dataset:
    """
    Function computing PCA on a block of a Dataset.
    Made for xarray.map_blocks.
    """

    if ds.dims[DIM_VARIANT] == 0: # then return template
        size_sample = ds.dims[DIM_SAMPLE]
        return xr.Dataset(
            data_vars={
                'vectors': ([DIM_SAMPLE, DIM_PC, DIM_WINDOW], np.zeros((size_sample, size_sample, 1), dtype=np.float64)),
                'vals': ([DIM_PC, DIM_WINDOW], np.zeros((size_sample, 1), dtype=np.float64)),
                'C': ([DIM_SAMPLE, DIM_SAMPLE, DIM_WINDOW], np.zeros((size_sample, size_sample, 1), dtype=np.float64)),
                'P': ([DIM_VARIANT, DIM_PC, DIM_WINDOW], np.zeros((window_size, size_sample, 1), dtype=np.float64)),
                'K': ([DIM_WINDOW], np.zeros((1), dtype=np.int8)),
            },
        )

    gl = ds.genotype_likelihood[:,:,:-1].values
    L = gl.reshape(gl.shape[0], -1)
    L = L.astype('float32')

    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
            f = shared.emMAF(L, iter=emMAF_iter, tole=emMAF_tole, t=emMAF_t)
            C, P, K = covariance.emPCA(L, f, e=emPCA_e, iter=emPCA_iter, tole=emPCA_tole, t=emPCA_t)
    vals, vectors = np.linalg.eig(C)

    pca_res = xr.Dataset(
        data_vars={
            'vectors': ([DIM_SAMPLE, DIM_PC, DIM_WINDOW], vectors[:,:,np.newaxis]),
            'vals': ([DIM_PC, DIM_WINDOW], vals[:,np.newaxis]),
            'C': ([DIM_SAMPLE, DIM_SAMPLE, DIM_WINDOW], C[:,:,np.newaxis]),
            'P': ([DIM_VARIANT, DIM_PC, DIM_WINDOW], P[:,:,np.newaxis]),
            'K': ([DIM_WINDOW], np.array(K, dtype=np.int8))
        },
    )
    return pca_res