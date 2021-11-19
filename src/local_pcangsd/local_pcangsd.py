import numpy as np
import pandas as pd
import xarray as xr
import sgkit as sg

DIM_VARIANT = "variants"
DIM_SAMPLE = "samples"
DIM_ALLELE = "alleles"
DIM_GENOTYPE = "genotypes"

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
    Wrapper arround sgkit.window_by_position
    Size either in bp or number of variants depending on type
    """
    if type == "position":
        return sg.window_by_position(ds, size=size)
    elif type == "variant":
        return sg.window_by_variant(ds, size=size)
    else:
        raise ValueError('window type accepted are only "position" or "variant"')



#====================================================

# function for calling pcangsd functions on dataset window