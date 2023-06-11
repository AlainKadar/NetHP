#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include "RandomNumbers.h"
#include "GraphCast.h"

class Integrator {
    private:

    bool added;
    int frame;
    int trial_i;
    int trial_j;


    public:

    RandomGenerator<float> rng_float;
    RandomGenerator<int> rng_int;
    
    dynamic_graph* N_ptr;
    
    bool initial_chunk;    
    
    bool final_chunk;    
    
    Integrator(dynamic_graph* N_ptr, int seed);

    ~Integrator();

    int rng(int upper_limit);
    
    float rng();

    // Degree constraint. Return true if violated
    bool const_deg(igraph_integer_t i, igraph_integer_t j);

    // Path constraint. Return true if violated
    bool const_path(igraph_integer_t i, igraph_integer_t j);
    
    // Planarity constraint. Return true if violated
    bool const_planar(igraph_integer_t i, igraph_integer_t j);

    // Pick a random vertex pair. If an edge exists, remove it. Otherwise,
    // add an edge between them, given that the degree constraint is not
    // violated
    // Return true if an edge toggle was proposed
    bool toggle_edge(igraph_integer_t i, igraph_integer_t j, bool print);

    // Return true if proposed move accepted
    bool transition();

    // Return graph energy
    float getEnergy(igraph_t *graph);

    void write_chunk(int timesteps, int dump_rate, int chunk_size,
            std::string write_dir);

};

#endif
