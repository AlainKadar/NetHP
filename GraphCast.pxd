# disutils: language = c++

from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "GraphCast.cpp":
    pass

cdef extern from "GraphCast.h" namespace "interface":
    cdef cppclass GraphCast:
        #void* N_ptr
        #igraph_t* G_ptr
        int num_verts
        float edge_cost
        float vert_cost
        float cluster_cost
        float temp
        int n
        float spikes
        int _len
        vector[int] rows
        vector[int] cols
        GraphCast() except +
        GraphCast(int, float, float, float, float) except +
        GraphCast(int, float, int, float, float, float, float) except +
        void Metropolis(int iterations, int chunk_size, int dump_rate, 
                        string write_dir)
        float getEnergy()
        void print_components()
        void getGraph()

