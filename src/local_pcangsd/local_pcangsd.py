#!python3
import gzip
import h5py
import numpy as np
import pandas as pd
from pcangsd_core import shared
from pcangsd_core import covariance
# import lostruct
import multiprocessing as mp
import os, contextlib

class local_pcangsd:
    """hdf5 class that stores the genotype likelihoods extracted for the beagle angsd format"""
    def __init__(self, file, pca_file=None):
        self.data = file
        try:
            with h5py.File(self.data) as data:
                self.n = len(data['samples'])
                self.chromosomes = [x for x in list(data.keys()) if x != 'samples']
        except:
            self.n = None
            self.chromosomes = None
        
        if pca_file is None:
            self.pca = None
        else:
            self.pca = pca_file
            if self.chromosomes is None:
                self.windows = None
            else:
                try:
                    with h5py.File(self.pca) as pca:
                        self.windows = {}
                        for chr in self.chromosomes:
                            self.windows[chr] = [(x[0], x[1]) for x in [y.split('-') for y in list(pca[chr].keys())]]
                except:
                    self.windows = None

    def populate(self, beagle_file, chunksize=10000, overwrite=False):
        """
        Given a beagle file from angsd, fill the data hdf5.
        """
        mode = 'w' if overwrite else 'x'
        try:
            with h5py.File(self.data, mode) as data:
                # read header
                with gzip.open(beagle_file, 'r') as fr:
                    ind_header = np.array(fr.readline().strip().split(b'\t')[3:])
                    _, idx = np.unique(ind_header, return_index=True)
                    samples = ind_header[np.sort(idx)]
                self.n = len(samples)
                data.create_dataset('samples', data=samples)

                # loop through beagle and add data
                kept_columns = [0, 1, 2] + [x for i, x in enumerate(range(3, (self.n*3)+3)) if (i+1) % 3 != 0]
                # this drops the 3rd genotype likelihood for each individual, same as new version of pcangsd
                for df in pd.read_csv(beagle_file, sep='\t', chunksize=chunksize, usecols=kept_columns):
                    df[['chr', 'pos']] = df.marker.str.rsplit('_', 1, expand=True)
                    add_to_hdf5(data, df)
                
                self.chromosomes = [x for x in list(data.keys()) if x != 'samples']
        except:
            print(
                "It seems the hdf5 file already exists.\n"
                "If you are unsure the populate method has completed successfully "
                "or if you want to overwrite, use the overwrite=True option"
                )

    def create_windows(self, win_size=10000, win_step=10000):
        """
        Create windows of snps over which local PCAs will be computed.
        Windows are based on number of SNPs for now.
        Indexing is python-like (0 based).
        """
        if self.chromosomes is None:
            print("Have you run the populate method?")
        with h5py.File(self.data, 'r') as data:
            self.windows = {}
            for chr in self.chromosomes:
                last_snp = data[chr]['pos'].shape[0]
                start = 0
                if last_snp <= win_size:
                    end = last_snp
                else:
                    end = win_size
                self.windows[chr] = [(start, end)]
                while end < last_snp:
                    start += win_step
                    end += win_step
                    self.windows[chr].append((start, end))
                last_start = start + win_step
                if last_start > last_snp:
                    self.windows[chr].append((last_start, last_snp))


    def compute_local_pca(self, max_cpu=1, pcangsd_cpu=1, overwrite=False):
        """
        Loop over all windows to compute the covariance matrices.
        This will erase previous results.
        The same K is used for all local PCAs to be compatible with lostruct.
        """
        n_cpu = min(mp.cpu_count(), max_cpu)

        if (overwrite==False) & (os.path.exists(self.pca)):
            print(f'File {self.pca} already exists, use option overwrite=True to replace.')
        else:
            pca = h5py.File(self.pca, 'w')
            data = h5py.File(self.data, 'r')

            for chr in self.chromosomes:
                # parallelize on windows for given chr
                chunk_dict = {f"{win[0]}-{win[1]}": data[chr]['gl'][win[0]:win[1]] for win in self.windows[chr]}
                pool = mp.Pool(n_cpu)
                result_pool = {win: pool.apply_async(emPCA_window, args=(L, 15, 0, 100, 1e5, 1)) for win, L in chunk_dict.items()}
                res_chr_dict = {win: r.get() for win, r in result_pool.items()}
                pool.close()
                pool.join()
                # save this chr to hdf5
                print(f'chromosome {chr} done')
                _ = pca.create_group(chr)
                for win in self.windows[chr]:
                    win_name = f"{win[0]}-{win[1]}"
                    _ = pca[chr].create_group(win_name)
                    _ = pca[chr][win_name].create_dataset('C', data=res_chr_dict[win_name][0])
                    _ = pca[chr][win_name].create_dataset('P', data=res_chr_dict[win_name][1])
                    _ = pca[chr][win_name].create_dataset('K', data=res_chr_dict[win_name][2])
                    _ = pca[chr][win_name].create_dataset('eigval', data=res_chr_dict[win_name][3])
                    _ = pca[chr][win_name].create_dataset('eigvec', data=res_chr_dict[win_name][4])
            
            data.close()
            pca.close()

    def print_dataset(self):
        with h5py.File(self.data, 'r') as data:
            groups = list(data.keys())
            print('\n'.join(groups))

    def print_pca(self):
        with h5py.File(self.pca, 'r') as pca:
            groups = list(pca.keys())
            for g in groups:
                print(f'{g}')
                for w in pca[g].keys():
                    print(f'\t{w}')

    def to_lostruct(self):
        """
        Convert the hdf5 pca results to an array understandable by lostruct.
        """
        with h5py.File(self.pca, 'r') as pca:
            results = []
            win_names = []
            for chr in self.chromosomes:
                for win in pca[chr].keys():
                    win_names.append(f"{chr}:{win}")
                    results.append((pca[chr][win]['C'][:], None, pca[chr][win]['eigval'][:], pca[chr][win]['eigvec'][:]))
        results = np.array(results, dtype=object)
        return results, win_names

#===================================================================
# Helper functions

def add_to_hdf5(hdf5, df):
    """
    Add a dataframe directly imported from beagle to the hdf5.
    One group per chromosome.
    """
    for chr, tmpdf in df.groupby('chr'):
        if chr not in hdf5.keys():
            _ = hdf5.create_group(chr)
        if 'pos' not in hdf5[chr].keys():
            # create dataset and add first data
            _ = hdf5[chr].create_dataset('pos', dtype='int', maxshape=(None,), chunks=True,
                data=tmpdf.pos.to_numpy(dtype='int'))
        else:
            # extend dataset
            x1 = hdf5[chr]['pos'].shape[0]
            x2 = x1 + tmpdf.shape[0]
            hdf5[chr]['pos'].shape = (x2,)
            hdf5[chr]['pos'][x1:x2] = tmpdf.pos.to_numpy(dtype='int')
        if 'gl' not in hdf5[chr].keys():
            # create dataset and add first data
            # pcangsd expects float32
            _ = hdf5[chr].create_dataset('gl', dtype='float32', maxshape=(None, hdf5['samples'].shape[0]*3), chunks=True,
                data=tmpdf.iloc[:,3:-2].to_numpy(dtype='float32')) # -2 for the two added columns at the end 'chr' and 'pos'
        else:
            # extend dataset
            x1, y = hdf5[chr]['gl'].shape
            x2 = x1 + tmpdf.shape[0]
            hdf5[chr]['gl'].shape = (x2, y)
            hdf5[chr]['gl'][x1:x2,:] = tmpdf.iloc[:,3:-2].to_numpy(dtype='float32')


def emPCA_window(L, global_K, e=0, iter=100, tole=1e-5, t=1):
    # for now do not care about min_MAF
    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
            f = shared.emMAF(L, iter=200, tole=1e-4, t=t)
            C, P, K = covariance.emPCA(L, f, e=e, iter=iter, tole=tole, t=t)
    vals, vectors = np.linalg.eig(C)
    vals = vals[:global_K].real.astype('float32')
    vectors = vectors[:, :global_K].T.real.astype('float32')
    return (C, P, global_K, vals, vectors)

