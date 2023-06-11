#include <igraph/igraph.h>
#include <string>

class dynamic_graph : public igraph_t {
    public:
    float edge_cost;
    float vert_cost;
    float cluster_cost;
    float temp;
    int num_verts;
    float current_energy;

    // Nullary constructor
    dynamic_graph();

    // Constructor - initialise an empty graph
    dynamic_graph(int num_verts, float edge_cost, float vert_cost,
            float cluster_cost, float temp);

    // Constructor - initialise a square lattice decorated with spikes
    dynamic_graph(int n, float spikes, int len, float edge_cost, float vert_cost,
            float cluster_cost, float temp);
    
    // Destructor
    ~dynamic_graph();

    // Print the component size distribution
    void print_components();

    void Metropolis(int timesteps, int dumprate, int chunk_size, 
            std::string write_dir, int seed);

    void write_report(std::string write_dir);
};



