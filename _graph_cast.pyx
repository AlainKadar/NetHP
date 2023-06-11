# distutils: language = c++

from GraphCast cimport GraphCast

import numpy as np
import networkx as nx
import scipy.sparse as sp
import os

from libcpp.string cimport string

cdef class PyCast:
    cdef GraphCast* c_cast

    # Nullary constructor
    def __cinit__(self):
        self.c_cast = new GraphCast()

    # Default constructor
    """
    def __cinit__(self, int num_verts, float edge_cost, float vert_cost,
                  float cluster_cost, float temp):

        self.c_cast = new GraphCast(num_verts, edge_cost, vert_cost,
                cluster_cost, temp)

        self.c_cast.num_verts = num_verts
    """

    # Spiky lattice constructor
    def __cinit__(self, float edge_cost, float vert_cost, float cluster_cost,
                  float temp, int n=0, float spikes=0, int _len=0,
                  stability=False):

        if stability:
            self.c_cast = new GraphCast(n, spikes, _len, edge_cost,
                                        vert_cost, cluster_cost, temp)
        else:
            self.c_cast = new GraphCast(n, edge_cost, vert_cost,
                                        cluster_cost, temp)

    #def __dealloc__(self):
    #    del self.c_cast
    
    def Metropolis(self, int timesteps, int chunk_size, int dump_rate,
                   string write_dir, int seed=1750):
        print("Running metropolis")
        if not os.path.isdir(write_dir):
            raise ValueError(write_dir, " is not a directory")

        self.c_cast.Metropolis(timesteps, chunk_size, dump_rate, write_dir,
                               seed)

    @property
    def getEnergy(self):
        return self.c_cast.getEnergy()

    def print_components(self):
        self.c_cast.print_components()

    @property
    def Graph(self):
        self.c_cast.getGraph()
        _len = len(self.c_cast.rows)
        matrix = sp.coo_matrix((np.ones((_len)),
                               (self.c_cast.rows, self.c_cast.cols)),
                               shape=(self.c_cast.num_verts, self.c_cast.num_verts))
        graph = nx.from_scipy_sparse_array(matrix)

        return graph

