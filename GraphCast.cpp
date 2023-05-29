#include "GraphCast.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>

#define MAX_DEGREE 4
#define MAX_PATH 3

/* PROGRAM FOR THE METRPOLIS SAMPLING OF EXPONENTIAL RANDOM GRAPHS
 *   
 *   Algorithm 1:
 *   Each Markov move is either the removal or addition of an
 *   edge
 *   The maximum degree per node is constrained
 *   See Algorithm 2 for edge swapping Markov moves
 *
 *   */

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

std::vector<int>& operator<<=(std::vector<int>& V, igraph_vector_int_t& I) {
    for (int i=0; i<igraph_vector_int_size(&I); i++) {
        V.push_back(VECTOR(I)[i]);
    }
    return V;
}

namespace interface {

    // Nullary constructor
    GraphCast::GraphCast () : N_ptr() {} 

    // Default constructor
    GraphCast::GraphCast (int num_verts, float edge_cost, float vert_cost,
            float cluster_cost, float temp) : N_ptr(num_verts, edge_cost, vert_cost, cluster_cost, temp)
    {
        //Network graph(num_verts, edge_cost, vert_cost, cluster_cost, temp);
        G_ptr = &N_ptr;
    }


    // Spiky lattice constructor
    GraphCast::GraphCast (int n, float spikes, int len, float edge_cost,
            float vert_cost, float cluster_cost, float temp) : N_ptr(n, spikes, len, edge_cost, vert_cost, cluster_cost,
                temp) {
        //Network graph(n, spikes, len, edge_cost, vert_cost, cluster_cost,
        //        temp);
        //N_ptr = graph;
        G_ptr = &N_ptr;
        num_verts = N_ptr.num_verts;
    }


    // Destructor
    GraphCast::~GraphCast () {
        igraph_destroy(&N_ptr);
    }

    void GraphCast::getGraph () {
        igraph_sparsemat_t adj_matrix;
        igraph_vector_t weights;
        igraph_vector_int_t i, j;
        igraph_vector_t x;
        
        igraph_vector_int_init(&i, igraph_ecount(&N_ptr));
        igraph_vector_int_init(&j, igraph_ecount(&N_ptr));
        igraph_vector_init(&weights, igraph_ecount(&N_ptr));
        igraph_vector_init(&x, igraph_ecount(&N_ptr));
        igraph_sparsemat_init(&adj_matrix, num_verts, num_verts,
                igraph_ecount(&N_ptr));

        igraph_get_adjacency_sparse(&N_ptr, &adj_matrix,
               IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_LOOPS_TWICE);
        igraph_sparsemat_getelements(&adj_matrix, &i, &j, &x);
        rows <<= i;
        cols <<= j;
        igraph_sparsemat_destroy(&adj_matrix);
        igraph_vector_int_destroy(&i);
        igraph_vector_int_destroy(&j);
        igraph_vector_destroy(&x);
        igraph_vector_destroy(&weights);
    }

    void GraphCast::Metropolis (int iterations, int chunk_size, int dump_rate,
            std::string write_dir) {
        this->N_ptr.Metropolis(iterations, chunk_size, dump_rate, write_dir);
    }
    
    float GraphCast::getEnergy () {
        return this->N_ptr.current_energy;
    }
    
    void GraphCast::print_components () {
        this->N_ptr.print_components();
    }
}

Integrator::Integrator(Network *N_ptr) {
    this->N_ptr = N_ptr;
    this->added = false;
    this->initial_chunk = true;
    this->final_chunk = false;
    this->N_ptr->current_energy = getEnergy(this->N_ptr);
}

Integrator::~Integrator()
{
    //delete this->N_ptr;
}
// Pick a random vertex pair. If an edge exists, remove it. Otherwise,
// add an edge between them, given that the degree constraint is not
// violated
// Return whether an edge was toggled or not
bool Integrator::toggle_edge(igraph_integer_t i, igraph_integer_t j) {
    igraph_integer_t eid;
    igraph_get_eid(this->N_ptr, &eid, i, j, false, false);
    if (eid==-1) { //i.e. no edge exists
        if (const_deg(i,j) or const_path(i,j) or const_planar(i,j)) { //check whether adding edge would violate constraints
            return false;
        };
        igraph_add_edge(this->N_ptr, i, j);
        this->added = true;
    } else { 
        igraph_delete_edges(this->N_ptr, igraph_ess_1(eid));
        this->added = false;
    }
    return true;
}

// Return true if proposed move accepted
bool Integrator::transition() {
    igraph_integer_t i = igraph_rng_get_integer(igraph_rng_default(),
            0, this->N_ptr->num_verts-1);
    igraph_integer_t j = igraph_rng_get_integer(igraph_rng_default(),
            0, this->N_ptr->num_verts-1);
    this->i=i;
    this->j=j;
    if (i==j) { return false; }
    
    if (toggle_edge(i,j)) {
        //std::cout << "Delta E is " << getEnergy(this)-this->current_energy << std::endl;
        if (getEnergy(this->N_ptr)>this->N_ptr->current_energy) {
            if (exp(-(getEnergy(this->N_ptr)-this->N_ptr->current_energy)/this->N_ptr->temp)<igraph_rng_get_unif01(igraph_rng_default())) {
                toggle_edge(i,j);
                //std::cout << "Reverse the move" << std::endl;
                return false;
            }
        }
        this->N_ptr->current_energy = getEnergy(this->N_ptr);
        return true;
    }
    return false;
}
// Returns true if vertex id is on the boundary of a square lattice
bool boundary(int i, int n) {
    if (i<n) { return true; }
    if (i % n == 0) { return true; }
    if ((i+1) % n == 0) { return true; }
    if (i>=(n*(n-1))) { return true; }
    return false;
}
// Constructor - default
Network::Network() {}

// Constructor - initialise an empty graph
Network::Network(int num_verts, float edge_cost, float vert_cost, float temp,
        float cluster_cost) {
    igraph_empty(this, num_verts, IGRAPH_UNDIRECTED);
    this->num_verts = num_verts;
    this->edge_cost = edge_cost;
    this->vert_cost = vert_cost;
    this->cluster_cost = cluster_cost;
    this->temp = temp;
}

    // Constructor - initialise a square lattice decorated with spikes
Network::Network(int n, float spikes, int len, float edge_cost, float vert_cost,
        float cluster_cost, float temp) {
    igraph_vector_int_t dims;
    igraph_vector_int_init(&dims, 2);
    VECTOR(dims)[0] = n;
    VECTOR(dims)[1] = n;
    igraph_square_lattice(this, &dims, 1, false, false, NULL);

    this->num_verts = n*n;
    this->edge_cost = edge_cost;
    this->vert_cost = vert_cost;
    this->cluster_cost = cluster_cost;
    this->temp = temp;

    int total_verts = this->num_verts;
    for (int i=0; i<n*n; i++) {
        if (boundary(i,n)) {
            if (igraph_rng_get_unif01(igraph_rng_default())<spikes) {
                igraph_add_vertices(this, len, NULL);
                int new_edge_source=i;
                for (int j=0; j<len; j++) {
                    igraph_add_edge(this, new_edge_source, total_verts);
                    total_verts++;
                    new_edge_source=total_verts-1;
                }
            }
        }
    }
    num_verts = total_verts;
    igraph_vector_int_destroy(&dims);
}

// Destructor
Network::~Network() {
    igraph_destroy(this);
}


// Print the component size distribution
void Network::print_components() {
    igraph_vector_int_t csize;
    igraph_integer_t i;
    igraph_vector_int_init(&csize, 0);
    igraph_connected_components(this, NULL, &csize, 0, IGRAPH_STRONG);

    for (i = 0; i < igraph_vector_int_size(&csize); i++) {
        printf("%i, ", int(VECTOR(csize)[i]));
        //std::cout << VECTOR(csize)[i] << ", ";
    }
    printf("\n");
    //std::cout << std::endl;
    igraph_vector_int_destroy(&csize);
}

// Return graph energy
float Integrator::getEnergy(igraph_t *graph) {
    float energy = 0;
    igraph_vector_int_t membership, csize, component_ids;
    igraph_integer_t i, j;
    igraph_vs_t vs;

    /*EDGE ENERGY CONTRIBUTION*/
    energy += (this->N_ptr->edge_cost + this->N_ptr->temp)*igraph_ecount(this->N_ptr);
    
   
    // get the connected components
    igraph_vector_int_init(&membership, 0);
    igraph_vector_int_init(&csize, 0);
    igraph_vector_int_init(&component_ids, 0);
    igraph_connected_components(this->N_ptr, &membership, &csize, 0, IGRAPH_STRONG);

    // iterate over the subgraphs
    for (i = 0; i < igraph_vector_int_size(&csize); i++) {
        /*COMPONENT SIZE ENERGY CONTRIBUTION*/
        energy += this->N_ptr->vert_cost*((VECTOR(csize)[i]-1)*
                (VECTOR(csize)[i]-1));

        /*CLUSTER ENERGY CONTRIBUTION*/
        // Requires subgraphs
        for (j = 0; j < this->N_ptr->num_verts; j++) {
            if (VECTOR(membership)[j]==i) {
                igraph_vector_int_push_back(&component_ids, j);
            }
        }
        igraph_vs_vector(&vs, &component_ids);
        igraph_t subgraph;
        igraph_induced_subgraph(this->N_ptr, &subgraph, vs,
                IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH);
        igraph_real_t res;
        igraph_transitivity_undirected(&subgraph, &res, IGRAPH_TRANSITIVITY_ZERO);
        energy += this->N_ptr->cluster_cost*res;
        igraph_destroy(&subgraph);
    }
    igraph_vector_int_destroy(&membership);
    igraph_vector_int_destroy(&csize);
    igraph_vector_int_destroy(&component_ids);

    return energy;
}

// Degree constraint. Return true if violated
bool Integrator::const_deg(igraph_integer_t i, igraph_integer_t j) {
    igraph_integer_t deg_i, deg_j;
    igraph_degree_1(this->N_ptr, &deg_i, i, IGRAPH_ALL, false);
    igraph_degree_1(this->N_ptr, &deg_j, j, IGRAPH_ALL, false);
    return (deg_i==MAX_DEGREE or deg_j==MAX_DEGREE);
}

// Path constraint. Return true if violated
bool Integrator::const_path(igraph_integer_t i, igraph_integer_t j) {
    igraph_matrix_t res;
    igraph_matrix_init(&res, 1, 1);
    igraph_distances(this->N_ptr, &res, igraph_vss_1(i), igraph_vss_1(j), IGRAPH_ALL);
    bool violated = (MATRIX(res,i,j)>=MAX_PATH);
    igraph_matrix_destroy(&res);
    return violated;
}

bool Integrator::const_planar(igraph_integer_t i, igraph_integer_t j) {
    // Create an empty BGL graph with the same number of vertices as the igraph.
    Graph graph(igraph_vcount(this->N_ptr));

    // Iterate over the edges in the igraph.
    igraph_es_t es;
    igraph_es_all(&es, IGRAPH_EDGEORDER_ID);
    igraph_eit_t eit;
    igraph_eit_create(this->N_ptr, es, &eit);

    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t from, to;
        igraph_edge(this->N_ptr, IGRAPH_EIT_GET(eit), &from, &to);

        // Add the same edge to the BGL graph.
        add_edge(int(from), int(to), graph);

        IGRAPH_EIT_NEXT(eit);
    }

    // Add proposed edge to BGL graph
    add_edge(int(i),int(j),graph);

    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);
    
    return !boyer_myrvold_planarity_test(graph);
}

void Integrator::write_chunk(int chunk_size, int dump_rate, int chunk_num,
        std::string write_dir) {
    int row_num;
    int num_nodes = igraph_vcount(this->N_ptr);
    int size = 1000;
    std::vector<int> s(size, -1);
    std::vector<int> t(size, -1);
    std::vector<int> start(size, -1);
    std::vector<int> end(size, -1);
    std::vector<int> energy(chunk_size, -1);
    std::string edge_fname = write_dir + "edges.csv";
    std::string node_fname = write_dir + "nodes.csv";
    std::string energy_fname = write_dir + "energy.csv";
    if (this->initial_chunk) {
        std::remove(edge_fname.c_str());
        std::remove(node_fname.c_str());
        std::remove(energy_fname.c_str());
        std::ofstream edges_file(edge_fname, std::ios::app);
        std::ofstream nodes_file(node_fname);
        std::ofstream energy_file(energy_fname, std::ios::app);
        this->frame = 1;
        int num_edges = igraph_ecount(this->N_ptr);
        row_num = num_edges + 1;
        edges_file << "Source;Target;timeset" << std::endl;
        edges_file.close(); 
        nodes_file << "ID;Start;End" << std::endl;
        nodes_file.close(); 
        energy_file << "Step,Energy" << std::endl;
        energy_file.close(); 
        igraph_integer_t from, to;
        for (int i=0; i<num_edges; i++) {
            igraph_edge(this->N_ptr, i, &from, &to);
            s[i]=from;
            t[i]=to;
            start[i]=0;
        }
    } else { row_num=0; }

    std::ofstream energy_file(energy_fname, std::ios::app);
    for (int i=0; i<chunk_size; i++) {
        if (i % dump_rate == 0) {
            energy_file << chunk_num*chunk_size + i << "," << this->N_ptr->current_energy << std::endl;
        }
        if (transition()) {
            //std::cout << "Added?: " << this->added << std::endl;
            if (row_num<size) {
                if (this->added) {
                    s[row_num] = this->i;
                    t[row_num] = this->j;
                    start[row_num] = this->frame;
                    row_num++;
                } else {
                    for (int r=row_num; r>0; r--) {
                        if ((s[r]==this->i and t[r]==this->j)
                         or (s[r]==this->j and t[r]==this->i)) {
                            if (start[r]==-1) {
                            } else {
                                if (end[r]==-1) {
                                    end[r]=this->frame;
                                    break;
                                } else {
                                    continue;
                                }
                            }
                        }
                    }
                }
                this->frame++;
            } else {
                printf("Row number exceeded size!\n");
            }
        }
    }
    std::ofstream edges_file(edge_fname, std::ios::app);
    std::ofstream nodes_file(node_fname);
    for (int i=0; i<row_num; i++) {
        if (end[i]==-1) { end[i] = this->frame; }
        edges_file << s[i] << ";" << t[i] << ";<[" << start[i] << ", " << end[i] << "]>" << std::endl;
    }
    edges_file.close(); 
    for (int i=0; i<num_nodes; i++) {
        nodes_file << i << ";" << 0 << ";" << this->frame << std::endl;
    }
}


void Network::Metropolis(int timesteps, int chunk_size, int dump_rate,
       std::string write_dir) {
    Integrator I(this);
    int chunks = timesteps/chunk_size;
    for (int i=0; i<chunks; i++) {
        if (i==chunks-1) { I.final_chunk = true; }
        I.write_chunk(chunk_size, dump_rate, i, write_dir);
        I.initial_chunk = false;
    }
}

int main() {
    interface::GraphCast GC(4, 0.5, 5, -0.1, 0.1, 0.1, 1);
    //interface::GraphCast GC(10, -0.1, 0.1, 0.1, 1);
    std::string dir = "test/";
    GC.Metropolis(100000, 10000, 100, dir);
    //Network graph(10, 0.5, 5, -0.1, 0.1, 0.1, 1);
}
