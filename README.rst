=====
NetHP
=====

Overview
========
Graph Markov Chain Monte Carlo methods using the Metropolis Algorithm. Algorithms implemented in C++ with optional Python API,
Used to obtain results in our paper `"Graph–Property Relationships for Complex Chiral Nanodendrimers" <https://doi.org/10.1021/acsnano.4c12964>`__. This is a method known as Exponential Random Graph Models (ERGMs) in statistics and sometimes referred to as Network Hamiltonians.

Installation guide
==================
**NetHP** can be built from source for *macOS* using a conda environment and `setup.py`. Alternative platforms/methods will require editing `setup.py` with the correct paths to the **igraph** and **eigen** libraries. Clone the repository, install the build dependencies then install **NetHP**:

.. code:: bash

   git clone https://github.com/AlainKadar/NetHP
   cd NetHP
   conda install igraph eigen scipy boost
   python setup.py install

Example script
==============
Here is a simple example script that first instantiates an ERGM and then runs a Metropolis algorithm for :code:`100000` timesteps. The full list of available parameters can be found in the source code.

.. code:: python

   from MCMC import _graph_cast

   root = b"/directory/to/write/"
   N = _graph_cast.PyCast(-4, 1.5, 0, 0, 50, 1)
   N.Metropolis(100000, 100000, 100, root+b'1/', 11)
