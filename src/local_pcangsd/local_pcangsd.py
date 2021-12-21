import numpy as np
import pandas as pd
import xarray as xr
import sgkit as sg
import dask.array as da
import dask
from sgkit.window import _get_chunked_windows, _sizes_to_start_offsets
from pcangsd_core import shared
from pcangsd_core import covariance
import os, contextlib
from typing import Optional
import shutil

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


def window(ds: xr.Dataset, type: str, size: int, min_variant_number: int=100) -> xr.Dataset:
    """
    Wrapper arround sgkit.window_by_[...].
    Size either in bp or number of variants depending on type.
    """

    if type not in {"position", "variant"}:
        raise ValueError(
            "Window type accepted are only 'position' or 'variant'."
        )

    if type == "position":
        ds = sg.window_by_position(ds, size=size)
    elif type == "variant":
        ds = sg.window_by_variant(ds, size=size)

    # only keep windows with enough variants
    kept = np.where((ds.window_stop.values - ds.window_start.values) > min_variant_number)[0]
    return ds.isel(windows=kept)

def _pcangsd_wrapper(
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
            # hide messages from Pcangsd
            f = shared.emMAF(L, iter=emMAF_iter, tole=emMAF_tole, t=emMAF_t)
            C, P, _ = covariance.emPCA(L, f, e=emPCA_e, iter=emPCA_iter, tole=emPCA_tole, t=emPCA_t)
    C = C.astype(np.float32)
    vals, vectors = np.linalg.eig(C)
    vals = vals.astype(np.float32)
    vectors = vectors.astype(np.float32)
    total_variance = np.sum(np.power(C, 2).flatten())
    return C, total_variance, vals[:k], vectors[:, :k].T


def _create_save_pca_result(res, window_index, window_contigs, window_starts, window_stops, sample_id, attrs, tmp_folder, restart):
    tmp_file = f"{tmp_folder}/contig{window_contigs[window_index]}:{window_starts[window_index]}-{window_stops[window_index]}.zarr"
    if (not restart) or (not os.path.exists(tmp_file)):
        window_ds = xr.Dataset(
                data_vars={
                    'sample_id': ([DIM_SAMPLE], sample_id),
                    'window_contig': ([DIM_WINDOW], np.array([window_contigs[window_index]], dtype=np.int64)),
                    'window_start': ([DIM_WINDOW], np.array([window_starts[window_index]], dtype=np.int64)),
                    'window_stop': ([DIM_WINDOW], np.array([window_stops[window_index]], dtype=np.int64)),
                    'C': ([DIM_WINDOW, f'{DIM_SAMPLE}_0', f'{DIM_SAMPLE}_1'], res[0][np.newaxis]),
                    'total_variance': ([DIM_WINDOW], np.array([res[1]])),
                    'vals': ([DIM_WINDOW, DIM_PC], res[2][np.newaxis]),
                    'vectors': ([DIM_WINDOW, DIM_PC, DIM_SAMPLE], res[3][np.newaxis]),
                },
                attrs=attrs,
            )
        window_ds.to_zarr(tmp_file, mode="w")
    return tmp_file


def pca_window(
    ds: xr.Dataset, 
    zarr_store: str, 
    output_chunsize: int=10000,
    k: Optional[int]=None,
    emMAF_iter: int=200, emMAF_tole: float=1e-4, emMAF_t: int=1, 
    emPCA_e: int=0, emPCA_iter: int=100, emPCA_tole: float=1e-5, emPCA_t: int=1,
    tmp_folder: str='/tmp/tmp_local_pcangsd',
    scheduler: str='processes', num_workers: Optional[int]=None,
    clean_tmp: bool=True,
    restart: bool=True,
) -> np.array:
    """
    Run PCAngsd on each window.
    Return results as np.array of objects (lostruct style).

    Code modified from sgkit.window_statistic

    Arguments
    -------------------
    restart: do not overwrite tmp zarr file and continue, restarting the analysis in case of error.
    """
    
    if 'window_start' not in ds or 'window_stop' not in ds:
        raise Exception("Variables 'window_start' and 'window_stop' not defined in the Dataset.")

    if k is None:
        k = ds.dims['samples']

    values = ds.genotype_likelihood
    window_contigs = ds.window_contig.values
    window_starts = ds.window_start.values
    window_stops = ds.window_stop.values

    values = da.asarray(values)

    window_lengths = window_stops - window_starts
    depth = np.max(window_lengths)

    # Dask will raise an error if the last chunk size is smaller than the depth
    # Workaround by rechunking to combine the last two chunks in first axis
    # See https://github.com/dask/dask/issues/6597
    if depth > values.chunks[0][-1]:
        chunk0 = values.chunks[0]
        new_chunk0 = tuple(list(chunk0[:-2]) + [chunk0[-2] + chunk0[-1]])
        values = values.rechunk({0: new_chunk0})

    chunks = values.chunks[0]

    rel_window_starts, windows_per_chunk = _get_chunked_windows(chunks, window_starts, window_stops)

    # Add depth for map_overlap
    rel_window_starts = rel_window_starts + depth
    rel_window_stops = rel_window_starts + window_lengths

    chunk_offsets = _sizes_to_start_offsets(windows_per_chunk)

    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)

    def blockwise_moving_stat(x: np.array, block_info=None):
        if block_info is None or len(block_info) == 0:
            return np.array([])
        chunk_number = block_info[0]["chunk-location"][0]
        chunk_offset_start = chunk_offsets[chunk_number]
        chunk_offset_stop = chunk_offsets[chunk_number + 1]
        chunk_window_starts = rel_window_starts[chunk_offset_start:chunk_offset_stop]
        chunk_window_stops = rel_window_stops[chunk_offset_start:chunk_offset_stop]
        zarr_list = []
        window_index = windows_per_chunk[:chunk_number].sum() # number of windows in all previous chunks
        for i,j in zip(chunk_window_starts, chunk_window_stops):
            res_ij = _pcangsd_wrapper(x[i:j], k=k,
                emMAF_iter=emMAF_iter, emMAF_tole=emMAF_tole, emMAF_t=emMAF_t,
                emPCA_e=emPCA_e, emPCA_iter=emPCA_iter, emPCA_tole=emPCA_tole, emPCA_t=emPCA_t,
            )
            tmp_file = _create_save_pca_result(res_ij, window_index, window_contigs, window_starts, 
                window_stops, ds.sample_id.values, ds.attrs, tmp_folder, restart)
            zarr_list.append(tmp_file)
            window_index += 1
        return np.array(zarr_list, dtype=object)

    map_depth = {0: depth}

    result = values.map_overlap(
        blockwise_moving_stat,
        depth=map_depth,
        dtype=object,
        chunks=(1,),
        drop_axis=(1, 2),
        boundary=0,
        trim=False,
    )

    zarr_list = result.compute(scheduler=scheduler, num_workers=num_workers)

    ds_pca = xr.open_mfdataset(zarr_list, combine='nested', concat_dim='windows', engine="zarr").chunk({'windows': output_chunsize})
    to_store = ds_pca.copy()
    for var in to_store.variables:
        to_store[var].encoding.clear()
    to_store.to_zarr(zarr_store, mode='w')

    if clean_tmp:
        shutil.rmtree(tmp_folder, ignore_errors=True)
    
    return zarr_store


def to_lostruct(ds_pca: xr.Dataset) -> np.array:
    vectors = ds_pca.vectors.values
    vals = ds_pca.vals.values
    C = ds_pca.C.values
    total_variance = ds_pca.total_variance.values
    
    results = [
        (c, t, val, vec) for c, t, val, vec in zip(C, total_variance, vals, vectors)
    ]
    return np.array(results, dtype=object)


def pcangsd_merged_windows(
        ds,
        windows_idx: np.array,
        k: Optional[int]=None,
        emMAF_iter: int=200, emMAF_tole: float=1e-4, emMAF_t: int=1, 
        emPCA_e: int=0, emPCA_iter: int=100, emPCA_tole: float=1e-5, emPCA_t: int=1
    ) -> tuple:
    """
    Compute PCAngsd on merged windows of interest.
    """

    if 'window_start' not in ds or 'window_stop' not in ds:
        raise Exception("Variables 'window_start' and 'window_stop' not defined in the Dataset.")

    if k is None:
        k = ds.dims['samples']

    window_start = ds.window_start.values[windows_idx]
    window_stop = ds.window_stop.values[windows_idx]
    variants_idx = np.concatenate([np.array(range(i, j)) for i, j in zip(window_start, window_stop)])

    with dask.config.set(**{'array.slicing.split_large_chunks': False}):
        result_pca = _pcangsd_wrapper(
            ds.genotype_likelihood.isel(variants=variants_idx).values, 
            k=k,
            emMAF_iter=emMAF_iter, emMAF_tole=emMAF_tole, emMAF_t=emMAF_t,
            emPCA_e=emPCA_e, emPCA_iter=emPCA_iter, emPCA_tole=emPCA_tole, emPCA_t=emPCA_t,
        )
    
    return result_pca
