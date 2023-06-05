#ifndef GRAPHCAST_H
#define GRAPHCAST_H

#include <vector>
#include <string>
#include <igraph/igraph.h>
#include <random>

class Network : public igraph_t {
    public:
    float edge_cost;
    float vert_cost;
    float cluster_cost;
    float temp;
    int num_verts;
    float current_energy;

    // Nullary constructor
    Network();

    // Constructor - initialise an empty graph
    Network(int num_verts, float edge_cost, float vert_cost,
            float cluster_cost, float temp);

    // Constructor - initialise a square lattice decorated with spikes
    Network(int n, float spikes, int len, float edge_cost, float vert_cost,
            float cluster_cost, float temp);
    
    // Destructor
    ~Network();

    // Print the component size distribution
    void print_components();

    void Metropolis(int timesteps, int dumprate, int chunk_size, 
            std::string write_dir);

    void write_report(std::string write_dir);
};


template<typename T>
class RandomGenerator;

template<>
class RandomGenerator<int> {
private:
    std::mt19937 gen;
    std::uniform_int_distribution<> dis;

public:
    // Default constructor
    RandomGenerator() : gen(std::random_device()()), dis(0, 100) {} 

    // Constructor with parameters
    RandomGenerator(int seed, int lower_bound, int upper_bound)
        : gen(seed), dis(lower_bound, upper_bound) {}

    // Generate a random number
    int getRandomNumber() {
        return dis(gen);
    }
};

template<>
class RandomGenerator<float> {
private:
    std::mt19937 gen;
    std::uniform_real_distribution<float> dis;

public:
    // Default constructor
    RandomGenerator() : gen(std::random_device()()), dis(0.0f, 1.0f) {} 

    // Constructor with parameters
    RandomGenerator(int seed, float lower_bound, float upper_bound)
        : gen(seed), dis(lower_bound, upper_bound) {}

    // Generate a random number
    float getRandomNumber() {
        return dis(gen);
    }
};

class Integrator {
    private:

    bool added;
    int frame;
    int trial_i;
    int trial_j;


    public:

    RandomGenerator<float> rng_float;
    RandomGenerator<int> rng_int;
    
    Network* N_ptr;
    
    bool initial_chunk;    
    
    bool final_chunk;    
    
    Integrator(Network* N_ptr);

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

namespace interface {
    class GraphCast {
        public:
            Network N_ptr;
            igraph_t* G_ptr;
            int num_verts;
            float edge_cost;
            float vert_cost;
            float cluster_cost;
            float temp;
            int n;
            float spikes;
            int _len;
            GraphCast();
            GraphCast(int num_verts, float edge_cost, float vert_cost,
                    float cluster_cost, float temp);
            GraphCast(int n, float spikes, int _len, float edge_cost,
                    float vert_cost, float cluster_cost, float temp);
            ~GraphCast();
            void getGraph();
            void Metropolis(int iterations, int chunk_size, int dump_rate, 
                    std::string write_dir);
            float getEnergy();
            void print_components();
            std::vector<int> rows;
            std::vector<int> cols;
    };
}

#endif
