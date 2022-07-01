import numpy as np
import pandas as pd
import xarray as xr
import sgkit as sg
import dask.array as da
import dask
from sgkit.window import _get_chunked_windows, _sizes_to_start_offsets
from pcangsd import shared
from pcangsd import covariance
import os
import contextlib
from typing import Optional
import shutil

DIM_VARIANT = "variants"
DIM_SAMPLE = "samples"
DIM_ALLELE = "alleles"
DIM_GENOTYPE = "genotypes"
DIM_PC = "PCs"
DIM_WINDOW = "windows"


def _create_dataset(df: pd.DataFrame, contigs: list) -> xr.Dataset:
    """Creates the genotype likelihood dataset from a chunk of dataframe.

    This is used to work on chunks in beagle_to_zarr.

    Args:
        df: pandas DataFrame (chunk of the beagle file)
        contigs: list of contig names

    Returns:
        xarray.Dataset: dataset of the chunk.
    """
    sample_id = np.array(df.columns[range(3, df.shape[1], 3)])

    variants = df.marker.str.rsplit("_", 1, expand=True)
    variants.columns = ["contig", "position"]

    flat_gl = df.iloc[:, 3:].to_numpy(dtype=np.float64)
    ncol = flat_gl.shape[1]
    GL = np.stack(
        [
            flat_gl[:, range(0, ncol, 3)],
            flat_gl[:, range(1, ncol, 3)],
            flat_gl[:, range(2, ncol, 3)],
        ],
        axis=-1,
    )

    variant_contig = np.array(
        [contigs.index(n) for n in variants.contig], dtype=np.int64
    )

    ds = xr.Dataset(
        data_vars={
            "variant_contig_name": ([DIM_VARIANT], variants.contig.to_numpy()),
            "variant_contig": ([DIM_VARIANT], variant_contig),
            "variant_position": (
                [DIM_VARIANT],
                variants.position.to_numpy(dtype=np.int64),
            ),
            "sample_id": ([DIM_SAMPLE], sample_id),
            "allele": (
                [DIM_VARIANT, DIM_ALLELE],
                df.iloc[:, [1, 2]].to_numpy(dtype=np.int8),
            ),
            "genotype_likelihood": ([DIM_VARIANT, DIM_SAMPLE, DIM_GENOTYPE], GL),
        },
        attrs=dict(contigs=contigs),
    )
    return ds


def beagle_to_zarr(input: str, store: str, chunksize: int = 10000) -> None:
    """Converts an ANGSD genotype likelihood dataset to a Zarr array on disk.

    Args:
        input: path to a genotype likelihood file in beagle format produced by ANGSD.
        store: output file on disk to store the dataset as a zarr file, ex: "output.zarr".
        chunksize: size of each chunk in the variant dimension.
    """

    # first pass to obtain contig names
    df_chunked = pd.read_csv(input, sep="\t", chunksize=chunksize)
    contigs = [
        df.marker.str.rsplit("_", 1, expand=True).iloc[:, 0].unique()
        for df in df_chunked
    ]
    contigs = list(np.unique(np.hstack(contigs)))

    df_chunked = pd.read_csv(input, sep="\t", chunksize=chunksize)
    for i, df in enumerate(df_chunked):
        ds = _create_dataset(df, contigs=contigs)
        if i == 0:
            ds.to_zarr(store, mode="w")
        else:
            ds.to_zarr(store, append_dim="variants")


def load_dataset(store: str, **kwargs) -> xr.Dataset:
    """Wrapper around xarray.open_zarr

    Args:
        store: path to zarr store
        **kwargs: keyword arguments passed to xarray.open_zarr

    Returns:
        xarray.Dataset: the opened dataset
    """
    ds = xr.open_zarr(store, **kwargs)
    return ds


def window(
    ds: xr.Dataset, type: str, size: int, min_variant_number: int = 100
) -> xr.Dataset:
    """Create windows on the dataset.

    Wrapper arround sgkit.window_by_[...].
    Size either in bp or number of variants depending on type.

    Args:
        ds: input dataset.
        type: 'position' or 'variant'.
            Create windows either using their position or using the variant number
        size: size of each window, either in bp or number of variants depending on type.
        min_variant_number: minimal number of variants to keep the window.
            Windows with less than min_variant_number variants are discarded.

    Returns:
        xarray.Dataset: ds input with appended windowing variables.

    Raises:
        ValueError: if type is not 'position' or 'variant'.
    """

    if type not in {"position", "variant"}:
        raise ValueError("Window type accepted are only 'position' or 'variant'.")

    if type == "position":
        ds = sg.window_by_position(ds, size=size)
    elif type == "variant":
        ds = sg.window_by_variant(ds, size=size)

    # only keep windows with enough variants
    kept = np.where(
        (ds.window_stop.values - ds.window_start.values) >= min_variant_number
    )[0]
    return ds.isel(windows=kept)


def _pcangsd_wrapper(
    gl: np.array,
    k: int,
    **kwargs,
) -> tuple[np.array, float, np.array, np.array]:
    """Wrapper around PCAngsd functions

    returning a tuple
    similar to what is returned by lostruct.eigen_windows.

    Args:
        gl: genotype likelihoods as formatted in the beagle_to_zarr output
            i.e. (variants, samples, genotypes=3)
        k: number of PCs to retain.
        **kwargs: Arbitrary keyword arguments.
            Allows to pass pcangsd specific arguments.

    Returns:
        tuple: (covariance matrix, total variance, eigen values, eigen vectors)
    """

    L = gl[:, :, :-1].reshape(gl.shape[0], -1)
    L = L.astype("float32")

    args_emMAF = dict(
        iter=kwargs["maf_iter"],
        tole=kwargs["maf_tole"],
        t=kwargs["threads"],
    )
    args_emPCA = dict(
        e=kwargs["n_eig"],
        iter=kwargs["iter"],
        tole=kwargs["tole"],
        t=kwargs["threads"],
    )

    with open(os.devnull, "w") as devnull:
        with contextlib.redirect_stdout(devnull):
            # hide messages from Pcangsd
            f = shared.emMAF(L, **args_emMAF)
            C, P, _ = covariance.emPCA(L, f, **args_emPCA)
    C = C.astype(np.float32)
    vals, vectors = np.linalg.eig(C)
    vals = vals.astype(np.float32)
    vectors = vectors.astype(np.float32)
    total_variance = np.sum(vals)
    # total_variance = np.sum(np.power(C, 2).flatten())
    return C, total_variance, vals[:k], vectors[:, :k].T


def _create_save_pca_result(
    res: tuple[np.array, float, np.array, np.array],
    window_index,
    window_contigs,
    window_starts,
    window_stops,
    sample_id,
    attrs,
    tmp_folder: str,
    restart: bool,
) -> str:
    """Saves the window pca result to file.

    Args:
        res: result of _pcangsd_wrapper.
        window_index: index of the current window.
        window_contigs: contig of the current window.
        window_starts: start variant indexes of current window.
        window_stops: stop variant indexes of current window.
        sample_id: passed from original dataset.
        attrs: passed from original dataset.
        tmp_folder: path of temporary folder to use.
        restart: should result be overwritten. Overwritten if False.

    Returns:
        str: path to the tmp window store.
    """
    tmp_file = (
        f"{tmp_folder}/contig{window_contigs[window_index]}:"
        f"{window_starts[window_index]}-{window_stops[window_index]}.zarr"
    )
    if (not restart) or (not os.path.exists(tmp_file)):
        window_ds = xr.Dataset(
            data_vars={
                "sample_id": ([DIM_SAMPLE], sample_id),
                "window_contig": (
                    [DIM_WINDOW],
                    np.array([window_contigs[window_index]], dtype=np.int64),
                ),
                "window_start": (
                    [DIM_WINDOW],
                    np.array([window_starts[window_index]], dtype=np.int64),
                ),
                "window_stop": (
                    [DIM_WINDOW],
                    np.array([window_stops[window_index]], dtype=np.int64),
                ),
                "C": (
                    [DIM_WINDOW, f"{DIM_SAMPLE}_0", f"{DIM_SAMPLE}_1"],
                    res[0][np.newaxis],
                ),
                "total_variance": ([DIM_WINDOW], np.array([res[1]])),
                "vals": ([DIM_WINDOW, DIM_PC], res[2][np.newaxis]),
                "vectors": ([DIM_WINDOW, DIM_PC, DIM_SAMPLE], res[3][np.newaxis]),
            },
            attrs=attrs,
        )
        window_ds.to_zarr(tmp_file, mode="w")
    return tmp_file


def pca_window(
    ds: xr.Dataset,
    zarr_store: str,
    output_chunksize: int = 10000,
    k: Optional[int] = None,
    tmp_folder: str = "/tmp/tmp_local_pcangsd",
    scheduler: str = "threads",
    num_workers: Optional[int] = None,
    clean_tmp: bool = True,
    restart: bool = False,
    maf_iter: int = 200,
    maf_tole: float = 1e-4,
    n_eig: int = 0,
    iter: int = 100,
    tole: float = 1e-5,
    pcangsd_threads: int = 1,
) -> np.array:
    """Run PCAngsd on each window.

    Args:
        ds: local_pcangsd dataset containing genotype likelihoods and windows.
        zarr_store: path to store the local_pcangsd results.
        output_chunksize: size of chunks for the xarray output.
        k: number of PCs to retain in the output.
            By default will keep all.
        tmp_folder: folder to use to store temporary results.
            Will be created if it does not exist. /tmp by default.
        scheduler: dask single-machine scheduler to use.
            'threads', 'processes' or 'synchronous'.
        num_workers: dask number of workers to use.
            Be careful to adapt pcangsd_threads and this argument accordingly.
        clean_tmp: should the temporary folder by emptied?
        restart: should the analysis be restarted using existing temporary files?
            If False, will overwrite existing tmp files.
        maf_iter: pcangsd maf_iter argument.
        maf_tole: pcangsd maf_tole argument.
        n_eig: pcangsd n_eig argument.
        iter: pcangsd iter argument.
        tole: pcangsd tole argument.
        pcangsd_threads: pcangsd threads argument.
            Be careful to adapt num_workers and this argument accordingly.

    Returns:
        str: Path to the created zarr_store containing each window pcangsd.

    Raises:
        Exception: if window variables does not exist in ds.
    """

    if "window_start" not in ds or "window_stop" not in ds:
        raise Exception(
            "Variables 'window_start' and 'window_stop' not defined in the Dataset."
        )

    if k is None:
        k = ds.dims["samples"]

    args_pcangsd = dict(
        maf_iter=maf_iter,
        maf_tole=maf_tole,
        n_eig=n_eig,
        iter=iter,
        tole=tole,
        threads=pcangsd_threads,
    )

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

    rel_window_starts, windows_per_chunk = _get_chunked_windows(
        chunks, window_starts, window_stops
    )

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
        window_index = windows_per_chunk[
            :chunk_number
        ].sum()  # number of windows in all previous chunks
        for i, j in zip(chunk_window_starts, chunk_window_stops):
            res_ij = _pcangsd_wrapper(
                x[i:j],
                k=k,
                **args_pcangsd,
            )
            tmp_file = _create_save_pca_result(
                res_ij,
                window_index,
                window_contigs,
                window_starts,
                window_stops,
                ds.sample_id.values,
                ds.attrs,
                tmp_folder,
                restart,
            )
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

    zarr_list = result.compute(
        scheduler=scheduler,
        num_workers=num_workers,
        threads_per_worker=1,
    )

    ds_pca = xr.open_mfdataset(
        zarr_list, combine="nested", concat_dim="windows", engine="zarr"
    ).chunk({"windows": output_chunksize})
    to_store = ds_pca.copy()
    for var in to_store.variables:
        to_store[var].encoding.clear()
    to_store.to_zarr(zarr_store, mode="w")

    if clean_tmp:
        shutil.rmtree(tmp_folder, ignore_errors=True)

    return zarr_store


def to_lostruct(ds_pca: xr.Dataset) -> np.array:
    """Converts the local_pcangsd result to lostruct format

    Args:
        ds_pca: local_pcangsd PCA result dataset.

    Returns:
        numpy.array: array in the lostruct format.
    """
    vectors = ds_pca.vectors.values
    vals = ds_pca.vals.values
    C = ds_pca.C.values
    total_variance = ds_pca.total_variance.values

    results = [
        (c, t, val, vec) for c, t, val, vec in zip(C, total_variance, vals, vectors)
    ]
    return np.array(results, dtype=object)


def pcangsd_merged_windows(
    ds: xr.Dataset,
    windows_idx: np.array,
    k: Optional[int] = None,
    maf_iter: int = 200,
    maf_tole: float = 1e-4,
    n_eig: int = 0,
    iter: int = 100,
    tole: float = 1e-5,
    pcangsd_threads: int = 1,
) -> tuple:
    """Compute PCAngsd on merged windows of interest.

    Args:
        ds: local_pcangsd dataset containing genotype likelihoods and windows.
        windows_idx: indexes of windows to merge.
        k: number of PCs to retain in the output.
            By default will keep all.
        maf_iter: pcangsd maf_iter argument.
        maf_tole: pcangsd maf_tole argument.
        n_eig: pcangsd n_eig argument.
        iter: pcangsd iter argument.
        tole: pcangsd tole argument.
        pcangsd_threads: pcangsd threads argument.

    Returns:
        tuple: (covariance matrix, total variance, eigen values, eigen vectors)

    Raises:
        Exception: if the dataset do not have windows variables.
    """

    if "window_start" not in ds or "window_stop" not in ds:
        raise Exception(
            "Variables 'window_start' and 'window_stop' not defined in the Dataset."
        )

    if k is None:
        k = ds.dims["samples"]

    args_pcangsd = dict(
        maf_iter=maf_iter,
        maf_tole=maf_tole,
        n_eig=n_eig,
        iter=iter,
        tole=tole,
        threads=pcangsd_threads,
    )

    window_start = ds.window_start.values[windows_idx]
    window_stop = ds.window_stop.values[windows_idx]
    variants_idx = np.concatenate(
        [np.array(range(i, j)) for i, j in zip(window_start, window_stop)]
    )

    with dask.config.set(**{"array.slicing.split_large_chunks": False}):
        result_pca = _pcangsd_wrapper(
            ds.genotype_likelihood.isel(variants=variants_idx).values,
            k=k,
            **args_pcangsd,
        )

    return result_pca
