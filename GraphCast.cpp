#include "GraphCast.h"
#include "Integrators.h"

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
        //dynamic_graph graph(num_verts, edge_cost, vert_cost, cluster_cost, temp);
        G_ptr = &N_ptr;
    }


    // Spiky lattice constructor
    GraphCast::GraphCast (int n, float spikes, int len, float edge_cost,
            float vert_cost, float cluster_cost, float temp) : N_ptr(n, spikes, len, edge_cost, vert_cost, cluster_cost,
                temp) {
        //dynamic_graph graph(n, spikes, len, edge_cost, vert_cost, cluster_cost,
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
            std::string write_dir, int seed) {
        this->N_ptr.Metropolis(iterations, chunk_size, dump_rate, write_dir,
                seed);
    }
    
    float GraphCast::getEnergy () {
        return this->N_ptr.current_energy;
    }
    
    void GraphCast::print_components () {
        this->N_ptr.print_components();
    }
}

Integrator::Integrator(dynamic_graph *N_ptr, int seed) {
    this->N_ptr = N_ptr;
    this->added = false;
    this->initial_chunk = true;
    this->final_chunk = false;
    this->N_ptr->current_energy = getEnergy(this->N_ptr);
    this->rng_float = RandomGenerator<float>(seed, 0.0f, 1.0f);
    this->rng_int = RandomGenerator<int>(seed, 0, this->N_ptr->num_verts-1);
}

Integrator::~Integrator()
{
    //delete this->N_ptr;
}
// Pick a random vertex pair. If an edge exists, remove it. Otherwise,
// add an edge between them, given that the degree constraint is not
// violated
// Return whether an edge was toggled or not
bool Integrator::toggle_edge(igraph_integer_t i, igraph_integer_t j, bool print) {
    igraph_integer_t eid;
    igraph_get_eid(this->N_ptr, &eid, i, j, false, false);
    if (eid==-1) { //i.e. no edge exists
        if (const_deg(i,j) or const_planar(i,j)) { //check whether adding edge would violate constraints
            return false;
        };
        igraph_add_edge(this->N_ptr, i, j);
        //if (print) {printf("Added\n"); }
        this->added = true;
    } else { 
        igraph_delete_edges(this->N_ptr, igraph_ess_1(eid));
        this->added = false;
        //if (print) {printf("Deleted\n"); }
    }
    return true;
}

// Return true if proposed move accepted
bool Integrator::transition() {
    //std::cout << "Before transition ecount: " << int(igraph_ecount(this->N_ptr)) << std::endl;
    /*
    igraph_integer_t i = igraph_rng_get_integer(igraph_rng_default(),
            0, this->N_ptr->num_verts-1);
    igraph_integer_t j = igraph_rng_get_integer(igraph_rng_default(),
            0, this->N_ptr->num_verts-1);
    */


    this->trial_i=this->rng_int.getRandomNumber();
    this->trial_j=this->rng_int.getRandomNumber();
    if (this->trial_i==this->trial_j) { return false; }
    
    if (toggle_edge(this->trial_i,this->trial_j,false)) {
        //std::cout << "Delta E is " << getEnergy(this)-this->current_energy << std::endl;
        if (getEnergy(this->N_ptr)>this->N_ptr->current_energy) {
            float num = this->rng_float.getRandomNumber();
            //std::cout << num << std::endl;
            if (exp(-(getEnergy(this->N_ptr)-this->N_ptr->current_energy)/this->N_ptr->temp)<num) {
                //std::cout << "Before reversal ecount: " << igraph_ecount(this->N_ptr) << std::endl;
                toggle_edge(this->trial_i,this->trial_j,true);
                //std::cout << "Reverse the move" << std::endl;
                //std::cout << "After reversal ecount: " << igraph_ecount(this->N_ptr) << std::endl;
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
dynamic_graph::dynamic_graph() {}

// Constructor - initialise an empty graph
dynamic_graph::dynamic_graph(int num_verts, float edge_cost, float vert_cost, float temp,
        float cluster_cost) {
    igraph_empty(this, num_verts, IGRAPH_UNDIRECTED);
    this->num_verts = num_verts;
    this->edge_cost = edge_cost;
    this->vert_cost = vert_cost;
    this->cluster_cost = cluster_cost;
    this->temp = temp;
}

    // Constructor - initialise a square lattice decorated with spikes
dynamic_graph::dynamic_graph(int n, float spikes, int len, float edge_cost, float vert_cost,
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
    this->num_verts = total_verts;
    igraph_vector_int_destroy(&dims);
}

// Destructor
dynamic_graph::~dynamic_graph() {
    igraph_destroy(this);
}


// Print the component size distribution
void dynamic_graph::print_components() {
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
    float cluster_contribution = 0;
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
        //energy += this->N_ptr->cluster_cost*res;
        cluster_contribution += this->N_ptr->cluster_cost*res;

        igraph_destroy(&subgraph);
    }
    energy += cluster_contribution/igraph_vector_int_size(&csize);
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
    int size = 100000;
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
        std::ofstream nodes_file(node_fname, std::ios::app);
        std::ofstream energy_file(energy_fname, std::ios::app);
        this->frame = 1;
        int num_edges = igraph_ecount(this->N_ptr);
        row_num = num_edges;
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
            if (row_num>=size) {
                s.push_back(-1);
                t.push_back(-1);
                start.push_back(-1);
                end.push_back(-1);
                size++;
            }
            if (this->added) {
                s[row_num] = this->trial_i;
                t[row_num] = this->trial_j;
                start[row_num] = this->frame;
                row_num++;
            } else {
                for (int r=row_num; r>=0; r--) {
                    if ((s[r]==this->trial_i and t[r]==this->trial_j)
                     or (s[r]==this->trial_j and t[r]==this->trial_i)) {
                        if (end[r]==-1) {
                            end[r]=this->frame;
                            break;
                        } else {
                            continue;
                        }
                    }
                }
            }
            this->frame++;
        }
    }
    std::cout << "EOF2 ecount " << int(igraph_ecount(this->N_ptr)) << std::endl;
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
    std::cout << "EOF ecount " << int(igraph_ecount(this->N_ptr)) << std::endl;
}

void dynamic_graph::write_report(std::string write_dir) {
    std::string report_fname = write_dir + "report.txt";
    std::cout << "Report ecount: " << igraph_ecount(this) << std::endl;
    //Degree sequence
    igraph_vector_int_t degrees;
    igraph_vs_t vs;
    igraph_vs_all(&vs);
    igraph_vector_int_init(&degrees, igraph_vcount(this));
    igraph_degree(this, &degrees, vs, IGRAPH_ALL, false);
    int degree_dist[7] = {0, 0, 0, 0, 0, 0};
    for (int i=0; i<igraph_vcount(this); i++) {
        degree_dist[VECTOR(degrees)[i]]++;
    }
    std::remove(report_fname.c_str());
    std::ofstream report_file(report_fname, std::ios::app);
    for (int i=0; i<6; i++) {
        report_file << i << ": " << degree_dist[i] << std::endl;
    }
}

void dynamic_graph::Metropolis(int timesteps, int chunk_size, int dump_rate,
       std::string write_dir, int seed) {
    Integrator I(this, seed);
    int chunks = timesteps/chunk_size;
    for (int i=0; i<chunks; i++) {
        if (i==chunks-1) { I.final_chunk = true; }
        I.write_chunk(chunk_size, dump_rate, i, write_dir);
        I.initial_chunk = false;
        std::cout << "Metropolis ecount: " << igraph_ecount(this) << std::endl;
    }
    write_report(write_dir);
    std::cout << "Final ecount: " << igraph_ecount(this) << std::endl;
}

int main() {
    //interface::GraphCast GC(3, 0, 0, -20, 0, 0, 2);
    interface::GraphCast GC(50, -1.5, 0.1, 0, 3);
    std::string dir = "assembly_test/";
    GC.Metropolis(100000, 100000, 100000, dir, 1750);
}
