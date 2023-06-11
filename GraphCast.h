#ifndef GRAPHCAST_H
#define GRAPHCAST_H

#include <vector>
#include "Networks.h"

namespace interface {
    class GraphCast {
        public:
            dynamic_graph N_ptr;
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
                    std::string write_dir, int seed);
            float getEnergy();
            void print_components();
            std::vector<int> rows;
            std::vector<int> cols;
    };
}

#endif
