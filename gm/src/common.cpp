//
// Created by lizhi on 10/31/20.
//
#include "../inc/common.h"
using namespace std;
bool compare_score(Score S1, Score S2){
    if(S1.score1 > S2.score1){
        return true;
    }else if(S1.score1 == S2.score1){
        if(S1.score2>S2.score2){
            return true;
        }else if(S1.score2 == S2.score2){
            if(S1.score3>S2.score3){
                return true;
            }else{
                return false;
            }
        }else{
            return false;
        }
    }else{
        return false;
    }
}
bool binary_search(unsigned int v,const unsigned int *array, unsigned int start, unsigned int end){
    for(unsigned int i=start;i<end;++i){
        if(array[i]==v){
            return true;
        }
    }
    return false;
}
bool compare_signature(unsigned int *sig1, unsigned int *sig2){
    for(unsigned int i=0;i<Signature_Properties;++i){
        if(sig2[i]<sig1[i]){
            return false;
        }
    }
    return true;
}
unsigned int get_score1(unsigned int *selected_nodes,unsigned int orders_len,unsigned int *neighbors,
                        unsigned int *offset,unsigned int v){
    unsigned int dout = 0;
    unsigned int din = 0;
    for(unsigned int idx = 0;idx < orders_len; ++idx){
        unsigned int node = selected_nodes[idx];
        if(binary_search(v,neighbors,offset[node],offset[node+1])){
            dout++;
        }
    }
    for(unsigned int idx = offset[v];idx<offset[v+1];idx++){
        unsigned int adj = neighbors[idx];
        if(binary_search(adj,selected_nodes,0,orders_len)){
            din++;
        }
    }
    return dout + din;
}
unsigned int get_score2(unsigned int *neighbors,unsigned int *selected_vertexes,unsigned int ordered_len,
                        unsigned int *offset, unsigned int v){
    unsigned int score = 0;
    std::set<unsigned int> non_visited_adjs;
    for(unsigned int idx=0;idx<ordered_len;++idx){
        unsigned int vertex = selected_vertexes[idx];
        for(unsigned int adj_idx = offset[vertex];adj_idx<offset[vertex+1];adj_idx++){
            unsigned int adj = neighbors[adj_idx];
            if(!binary_search(adj,selected_vertexes,0,ordered_len)){
                non_visited_adjs.insert(adj);
            }
        }
    }
    for(auto adj:non_visited_adjs){
        if(binary_search(adj,neighbors,offset[v],offset[v+1])){
            score++;
        }
    }
    return score;
}

unsigned int get_score3(const unsigned int *selected_vertexes,unsigned int sequence_len,
                        unsigned int *neighbors,unsigned int *offset, unsigned int v,unsigned int V){
    unsigned int score = 0;
    std::set<unsigned int> reachable_adjs;
    for(unsigned int idx=0;idx<sequence_len;++idx){
        unsigned int vertex = selected_vertexes[idx];
        reachable_adjs.insert(vertex);
        for(unsigned int adj_idx=offset[vertex];adj_idx<offset[vertex+1];adj_idx++){
            reachable_adjs.insert(neighbors[adj_idx]);
        }
    }
    std::set<unsigned int> non_reachable_vertex;
    for(unsigned int vertex=0;vertex<V;vertex++){
        if(reachable_adjs.find(vertex)==reachable_adjs.end()){
            non_reachable_vertex.insert(vertex);
        }
    }
    for(auto vertex: non_reachable_vertex){
        if(binary_search(vertex,neighbors,offset[v],offset[v+1])){
            score++;
        }
    }
    return score;
}
Graph::Graph(unsigned int mode,std::string input_file){
    std::ifstream infile;
    infile.open(input_file);
    if(!infile){
       cout<<"reading graph file error "<<input_file<<endl;
       exit(-1);
    }
    std::string line;
    std::vector<std::string> lines;
    const std::string delimter = "\t";
	double load_start = omp_get_wtime();
	unsigned int line_index = 0;
	unsigned int query_nodes;
    while(getline(infile,line)){
        if(line_index == 0){
            query_nodes = stoi(line);
            line_index++;
            continue;
        }
        line_index++;
        auto pos = line.find(delimter);
        if(pos == std::string::npos){
            continue;
        }
        lines.push_back(line);
    }
    infile.close();
    E = line_index - 1;
    V = query_nodes;
	double buildCSRs = omp_get_wtime();
    signatures = (unsigned int *)malloc(sizeof(unsigned int)*V*Signature_Properties);
    order_sequence = (unsigned int *)malloc(sizeof(unsigned int)*V);
    for(int i=0;i<V*Signature_Properties;++i){
        signatures[i] = 0;
    }
    vertexes = (unsigned int *)malloc(sizeof(unsigned int)*V);
    auto *neighbor_count = (unsigned int *)malloc(sizeof(unsigned int)*V);
    auto *neighbor_count2 = (unsigned int *)malloc(sizeof(unsigned int)*V);
    auto *temp_neighbors_len = (unsigned int *)malloc(sizeof(unsigned int)*V);
    auto *reverse_neighbors_len = (unsigned int *)malloc(sizeof(unsigned int)*V);
    for(unsigned int v=0;v<V;++v){
        vertexes[v] = v;
        temp_neighbors_len[v] = 0;
        neighbor_count[v] = 0;
        neighbor_count2[v] = 0;
        order_sequence[v] = 0;
        reverse_neighbors_len[v] = 0;
    }
    double build_signature_cost = 0.0f;
    for(auto it=lines.begin();it<lines.end();++it) {
        line = *it;
        auto pos = line.find(delimter);
        int start = stoi(line.substr(0, pos));
        int dst = stoi(line.substr(pos + 1, line.size() - pos - 1));
        temp_neighbors_len[start]++;
        reverse_neighbors_len[dst]++;
		double start_t = omp_get_wtime();
        signatures[start*Signature_Properties+Out_degree_offset]++;
        signatures[dst*Signature_Properties+In_degree_offset]++;
		double end_t = omp_get_wtime();
		build_signature_cost += (end_t - start_t);
    }
    neighbors = (unsigned int *)malloc(sizeof(unsigned int)*E);
    reverse_neighbors = (unsigned int *)malloc(sizeof(unsigned int)*E);
    offset = (unsigned int *)malloc(sizeof(unsigned int)*(V+1));
    offset2 = (unsigned int *)malloc(sizeof(unsigned int)*(V+1));
    offset[0] = 0;
    offset2[0] = 0;
    for(unsigned int v=0;v<V;++v){
        offset[v+1] = temp_neighbors_len[v]+offset[v];
        offset2[v+1] = reverse_neighbors_len[v] + offset2[v];
    }
    for(auto it=lines.begin();it<lines.end();++it) {
        line = *it;
        auto pos = line.find(delimter);
        int start = stoi(line.substr(0, pos));
        int dst = stoi(line.substr(pos + 1, line.size() - pos - 1));
        neighbors[offset[start]+neighbor_count[start]] = dst;
        reverse_neighbors[offset2[dst]+neighbor_count2[dst]] = start;
        neighbor_count[start]++;
        neighbor_count2[dst]++;
    }
    for(unsigned int v=0;v<V;++v){
        vector<unsigned int> temp_vec(neighbors+offset[v],neighbors+offset[v+1]);
        sort(temp_vec.begin(),temp_vec.end());
        for(int i=0;i<temp_vec.size();++i){
           neighbors[offset[v]+i] = temp_vec[i];
        }
    }
    for(unsigned int v=0;v<V;++v){
        vector<unsigned int> temp_vec(reverse_neighbors+offset2[v],reverse_neighbors+offset2[v+1]);
        sort(temp_vec.begin(),temp_vec.end());
        for(int i=0;i<temp_vec.size();++i){
           reverse_neighbors[offset2[v]+i] = temp_vec[i];
        }
    }
    double buildCSRe = omp_get_wtime();
    free(temp_neighbors_len);
    free(neighbor_count);
    double start = omp_get_wtime();
    if(!mode){
        sort_search_order();
        //find_leaders();
    }
}
void Graph::detect_one_hop_loop() const{
    for(unsigned int v=0;v<V;++v){
        for(unsigned int adj_idx=offset[v];adj_idx<offset[v+1];adj_idx++){
            unsigned int adj = neighbors[adj_idx];
            if(binary_search(v,neighbors,offset[adj],offset[adj+1])){
                signatures[v*Signature_Properties+one_hop_loop_counting_offset]++;
            }
        }
    }
}
void Graph::find_leaders(){
    unsigned int max_degree = 0;
    for(unsigned int i=0;i<V;++i){
        if(max_degree<(offset[i+1] - offset[i] + offset2[i+1] - offset2[i])){
            max_degree = offset[i+1] - offset[i] + offset2[i+1] - offset2[i];
        }
    }
    for(unsigned int i=0;i<V;++i){
        unsigned int v = order_sequence[i];
        if(max_degree == offset[v+1] - offset[v] + offset2[v+1] - offset2[v]){
            leader_vertexes.push_back(i);
        }
    }
    std::cout<<leader_vertexes[0]<<" "<<order_sequence[0]<<endl;
}
void Graph::detect_triangle() const {
    for(unsigned int v=0;v<V;++v){
        for(unsigned int adj_idx=offset[v];adj_idx<offset[v+1];adj_idx++){
            unsigned int adj = neighbors[adj_idx];
            for(unsigned int next_adj_idx=offset[adj];next_adj_idx<offset[adj+1];next_adj_idx++){
                unsigned int next_adj = neighbors[next_adj_idx];
                if(binary_search(v,neighbors,offset[next_adj],offset[next_adj+1])){
                    signatures[v*Signature_Properties+two_hop_loop_counting_offset]++;
                }
            }
        }
    }
}
void Graph::sort_search_order(){
    unsigned int max_out_degree = 0;
    unsigned int idx;
    parents_offset = (unsigned int *)malloc(sizeof(unsigned int)*(V+1));
    succs_offset = (unsigned int *)malloc(sizeof(unsigned int)*(V+1));
    succs_offset[0] = 0;
    parents_offset[0] = 0;
    auto * parent_count = (unsigned int *)malloc(sizeof(unsigned int)*V);
    auto * temp_parents_count = (unsigned int *)malloc(sizeof(unsigned int)*V);
    auto * succs_count = (unsigned int *)malloc(sizeof(unsigned int)*V);
    auto * temp_succs_count = (unsigned int *)malloc(sizeof(unsigned int)*V);
    for(unsigned int v=0;v<V;++v){
        parents_offset[v+1] = 0;
        succs_offset[v+1] = 0;
        parent_count[v] = 0;
        temp_parents_count[v] = 0;
        succs_count[v] = 0;
        temp_succs_count[v] = 0;
        if(max_out_degree < signatures[v*Signature_Properties+Out_degree_offset]){
            max_out_degree = signatures[v*Signature_Properties+Out_degree_offset];
            idx = v;
        }
    }
    order_sequence[0] = idx;
    unsigned int inserted_vertexes = 1;
    while(inserted_vertexes < V){
        Score max_score(0,0,0);
        for(unsigned int v = 0;v<V;++v){
            if(binary_search(v,order_sequence,0,inserted_vertexes)){
                continue;
            }
            unsigned int score1 = get_score1(order_sequence,inserted_vertexes,neighbors,offset,v);
            unsigned int score2 = get_score2(neighbors,order_sequence,inserted_vertexes,offset,v);
            unsigned int score3 = get_score3(order_sequence,inserted_vertexes,neighbors,offset,v,V);
            Score temp_score(score1,score2,score3);
            if(compare_score(temp_score,max_score)){
                max_score = temp_score;
                idx = v;
            }
        }
        order_sequence[inserted_vertexes++] = idx;
    }
    for(unsigned int v=0;v<V-1;++v){
        unsigned int current_vertex = order_sequence[v];
        for(unsigned int next_id=v+1;next_id<V;next_id++){
            unsigned int vertex = order_sequence[next_id];
            if(binary_search(vertex,neighbors,offset[current_vertex],offset[current_vertex+1])){
                parent_count[vertex]++;
            }
        }
    }
    for(unsigned int v=0;v<V;++v){
        parents_offset[v+1] = parents_offset[v] + parent_count[v];
    }
    parents = (unsigned int *)malloc(sizeof(unsigned int)*parents_offset[V]);
    for(unsigned int v=0;v<V-1;++v){
        unsigned int current_vertex = order_sequence[v];
        for(unsigned int next_id=v+1;next_id<V;next_id++){
            unsigned int vertex = order_sequence[next_id];
            if(binary_search(vertex,neighbors,offset[current_vertex],offset[current_vertex+1])){
                parents[parents_offset[vertex]+temp_parents_count[vertex]++] = v;
            }
        }
    }
    for(unsigned int v=0;v<V;++v){
        if(v==0){
            succs_count[v] = 0;
            continue;
        }
        unsigned int curr_vertex = order_sequence[v];
        for(idx=0;idx<v;++idx){
            unsigned int pre_vertex = order_sequence[idx];
            if(binary_search(pre_vertex,neighbors,offset[curr_vertex],offset[curr_vertex+1])){
                succs_count[curr_vertex]++;
            }
        }
    }
    for(unsigned int v=0;v<V;++v){
        succs_offset[v+1] = succs_offset[v] + succs_count[v];
    }
    succs = (unsigned int *)malloc(sizeof(unsigned int)*succs_offset[V]);
    for(unsigned int v=0;v<V;++v){
        if(v==0){
            continue;
        }
        unsigned int curr_vertex = order_sequence[v];
        for(idx=0;idx<v;++idx){
            unsigned int pre_vertex = order_sequence[idx];
            if(binary_search(pre_vertex,neighbors,offset[curr_vertex],offset[curr_vertex+1])){
                succs[succs_offset[curr_vertex]+temp_succs_count[curr_vertex]++] = idx;
            }
        }
    }
    free(parent_count);
    free(succs_count);
    free(temp_succs_count);
    free(temp_parents_count);
}
