#!python3

try:
    import argparse
    import local_pcangsd as lp
    import os
    import sys
except ImportError:
    raise ImportError("Something went wrong in importing modules")

parser = argparse.ArgumentParser(prog="local_pcangsd", description="Runs local_pcangsd")
parser.add_argument(
    "zarr",
    metavar="ZARR_GL",
    type=str,
    help="Zarr file for storing genotype likelihoods. Can be an input or intermediate output.",
)
parser.add_argument(
    "output", metavar="PCA_OUTPUT", type=str, help="Zarr output file of the PCA results.",
)
parser.add_argument(
    "-b",
    "--beagle",
    type=str,
    help="Genotype likelihood file in beagle format produced by ANGSD",
)
parser.add_argument(
    "--overwrite-zarr", action="store_true", help="Overwrite already existing zarr store.",
    dest="overwrite_zarr",
)
parser.add_argument(
    "-c", "--chunksize", type=int, help="Chunksize for the zarr store.",
    default=1000000,
)
parser.add_argument(
    '-wt', "--window-type", type=str,
    help="Type of windowing performed. Either 'position' or 'variant'",
    default="variant",
    choices=['position', 'variant'],
    dest='window_type',
)
parser.add_argument(
    '-ws', "--window-size", type=int,
    help="Number of bases or variants per window, depending on the --window-type argument.",
    default=10000,
    dest='window_size',
)
parser.add_argument(
    "--min-var-number", type=int,
    help="Minimum number of variants per window. Windows will less variants will be discarded.",
    default=500,
    dest='min_var_number',
)
parser.add_argument(
    "--tmp-folder", metavar="FOLDER", type=str,
    help="Temporary folder to store PCA results",
    default='/tmp/tmp_local_pcangsd',
    dest='tmp_folder',
)
parser.add_argument(
    "-k", type=int,
    help="Number of PCs to retain in the results",
    default=10,
)


def main():
    args = parser.parse_args()

    if os.path.exists(args.zarr):
        if args.overwrite_zarr:
            lp.beagle_to_zarr(input=args.beagle, store=args.zarr, chunksize=args.chunksize)
        else:
            print("zarr file already exists, will use it as input.")
            print("Use option --overwrite-zarr to force its re-creation.")
    else:
        if args.beagle is None:
            print("Please provide a genotype likelihood file with the --beagle option.")
            sys.exit(1)
        else:
            lp.beagle_to_zarr(input=args.beagle, store=args.zarr, chunksize=args.chunksize)

    ds = lp.load_dataset(args.zarr, chunks=args.chunksize)

    ds = lp.window(
        ds,
        type=args.window_type,
        size=args.window_size,
        min_variant_number=args.min_var_number
    )

    pca_zarr_store = lp.pca_window(
        ds,
        zarr_store=args.output,
        tmp_folder=args.tmp_folder,
        k=args.k,
    )

    print(f"local PCA results have been store in {pca_zarr_store}")


if __name__ == "__main__":
    main()
