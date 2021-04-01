//
// Created by lizhi on 10/31/20.
//

#ifndef GRAPH_MATCH_GPU_UGLY_COMMON_H
#define GRAPH_MATCH_GPU_UGLY_COMMON_H
#define Signature_Properties 2
#define In_degree_offset 0
#define Out_degree_offset 1
#define one_hop_loop_counting_offset 2
#define two_hop_loop_counting_offset 3
#define BLK_NUMS 48
#define BLK_DIM 512
#define WARPS_EACH_BLK (BLK_DIM/32)
#define WORK_UNITS (BLK_NUMS*WARPS_EACH_BLK)
#define GPU_TABLE_LIMIT 6000000000
#define HelperSize 15000
#define QUERY_NODES 12
#define MAX_NE 100
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <map>
#include <utility>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <unordered_map>
#include <stack>
#include <deque>
#include <random>
#include <cuda.h>
#include "omp.h"
class Graph{
public:
    unsigned int V;
    unsigned int E;
    unsigned int * vertexes;
    unsigned int * neighbors;
    unsigned int * reverse_neighbors;
    unsigned int * offset;
    unsigned int * offset2;
    unsigned int * parents_offset;
    unsigned int * parents;
    unsigned int * order_sequence;
    unsigned int * signatures;
    unsigned int * succs;
    unsigned int * succs_offset;
    unsigned int * c_sets;
    unsigned int * c_sets_offset;
    unsigned int * pattern_sig;
    std::vector<unsigned int> leader_vertexes;
    void detect_one_hop_loop() const;
    void detect_triangle() const;
    void sort_search_order();
    void find_leaders();
    Graph(unsigned int mode,std::string input_file);
};
class Score{
public:
    unsigned int score1;
    unsigned int score2;
    unsigned int score3;
    Score(unsigned int a, unsigned int b, unsigned int c){
        score1 = a;
        score2 = b;
        score3 = c;
    }
};
bool compare_score(Score S1, Score S2);
bool binary_search(unsigned int v,const unsigned int *array, unsigned int start, unsigned int end);
bool compare_signature(unsigned int *sig1, unsigned int *sig2);
unsigned int get_score1(unsigned int *selected_nodes,unsigned int orders_len,unsigned int *neighbors,
                        unsigned int *offset,unsigned int v);
unsigned int get_score2(unsigned int *neighbors,unsigned int *selected_vertexes,unsigned int ordered_len,
                        unsigned int *offset, unsigned int v);
unsigned int get_score3(const unsigned int *selected_vertexes,unsigned int sequence_len,
                        unsigned int *neighbors,unsigned int *offset, unsigned int v,unsigned int V);
#endif //GRAPH_MATCH_GPU_UGLY_COMMON_H
