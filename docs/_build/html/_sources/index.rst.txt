.. local_pcangsd documentation master file, created by
   sphinx-quickstart on Thu Jun 23 11:42:33 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for local_pcangsd
===============================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Installation
============

The easiest way to install dependencies and ``local_pcangsd`` is with a conda environment defined in the
`yaml recipe <https://github.com/alxsimon/local_pcangsd/blob/main/conda_env.yaml>`_

.. code-block:: console

   conda env create -f conda_env.yaml
   conda activate local_pcangsd
   pip install git+https://github.com/alxsimon/local_pcangsd.git

Quickstart
==========

#. Transform your genotype likelihood file with ``local_pcangsd.beagle_to_zarr()``
#. Open your converted dataset using ``local_pcangsd.load_dataset()``
#. Create windows with ``local_pcangsd.window()``
#. Run PCAngsd on each window with ``local_pcangsd.pca_window()``
#. Open the PCA dataset using ``local_pcangsd.load_dataset()``
#. Convert to the python lostruct format using ``local_pcangsd.to_lostruct()``
#. Use the data any way you want.

Please see this `example <https://github.com/alxsimon/local_pcangsd/blob/main/example.ipynb>`_
for more details.

API
===

Details on functions and their arguments.

.. automodule:: local_pcangsd
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
