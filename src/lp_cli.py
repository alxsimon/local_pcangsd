#!python3

try:
    import argparse
    import local_pcangsd as lp
    import os
except ImportError:
    raise ImportError("Something went wrong in importing modules")


parser = argparse.ArgumentParser(prog="local_pcangsd", description="Runs local_pcangsd")
parser.add_argument(
    "-b",
    "--beagle",
    metavar="FILE",
    help="Genotype likelihood file in beagle format produced by ANGSD",
)
parser.add_argument(
    "-z",
    "--zarr",
    metavar="FILE",
    help="Zarr file for storing genotype likelihoods. Can be an input or intermediate output.",
)
parser.add_argument(
    "--overwrite-zarr", action="store_true", help="Overwrite already existing zarr store.",
    dest="overwrite_zarr",
)
parser.add_argument(
    "-o", "--output", metavar="FILE", help="Zarr output file of the PCA results."
)
parser.add_argument(
    "-c", "--chunksize", metavar="FILE", help="Chunksize for the zarr store."
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
        lp.beagle_to_zarr(input=args.beagle, store=args.zarr, chunksize=args.chunksize)


if __name__ == "__main__":
    main()
