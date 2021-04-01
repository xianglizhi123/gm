#include <bits/stdc++.h>
#include "./inc/common.h"
#include <time.h>
using namespace std;
string getFileBaseName(string s) {
    char sep = '/';
    size_t i = s.rfind(sep, s.length());
    if (i != string::npos) {
        return(s.substr(i+1, s.length() - i));
    }
    return("");
}
void shuffle_array(unsigned int arr[], unsigned long long int n)
{
    // To obtain a time-based seed
    srand (time(NULL));
    unsigned seed = rand() % 100 + 1;;
    // Shuffling our array
    shuffle(arr, arr + n,
            default_random_engine(seed));
}

inline void chkerr(cudaError_t code)
{
    if (code != cudaSuccess)
    {
        cout<<cudaGetErrorString(code)<<endl;
        exit(-1);
    }
}
__device__ bool compare_signature_gpu(unsigned int *sig1, unsigned int *sig2){
    for(unsigned int i=0;i<Signature_Properties;++i){
        if(sig2[i]<sig1[i]){
            return false;
        }
    }
    return true;
}
__device__ bool basic_search(unsigned int node, unsigned int *buffer, unsigned int len){
    for(unsigned int idx = 0;idx<len;idx++){
        if(node == buffer[idx]){
            return true;
        }
    }
    return false;
}
__device__ void memory_copy(unsigned int *from,unsigned int *to,unsigned int lane_id,unsigned int copy_len){
    for(unsigned int idx=lane_id;idx<copy_len;idx+=32){
        to[idx] = from[idx];
    }
}
__device__ bool compare_sig(unsigned int sig1,unsigned int sig2){
    /*for(int i=0;i<1;++i){
        unsigned int a = sig1[i];
        unsigned int b = sig2[i];
        if(a&b!=a){
            return false;
        }
     }
    return true;*/
    return ((sig1&sig2) == sig1);
}

__global__ void initialize_searching(unsigned int *qSignatures,unsigned int *dSignatures,
                                     unsigned int *result_table,unsigned int *orderSequence,unsigned int U,
                                     unsigned long long int *lengths){
    __shared__ unsigned int sharedSignatures[Signature_Properties];
    unsigned int v = orderSequence[0];
    for(unsigned int i=threadIdx.x;i<Signature_Properties;i+=1024){
        sharedSignatures[i] = qSignatures[v*Signature_Properties+i];
    }
    __syncthreads();
    unsigned int globalIndex = threadIdx.x + blockIdx.x*1024;
    unsigned int totalThreads = 68*1024;
    //compare_signature_gpu(sharedSignatures,&dSignatures[u*Signature_Properties])&&
    for(unsigned int u=globalIndex;u<U;u+=totalThreads){
        if(compare_signature_gpu(sharedSignatures,&dSignatures[u*Signature_Properties])){
            //if(compare_signature_gpu(sharedSignatures,&dSignatures[u*Signature_Properties])){
            unsigned int index = atomicAdd(&lengths[0],1);
            result_table[index] = u;
        }
    }
}
__device__ void initialize_intersection(unsigned int *pre_copy,unsigned int first_index,char join_type,unsigned int * intersection,
                                        unsigned int *helper_buffer,char *intersection_bits,unsigned int *offset,unsigned int *offset2,
                                        unsigned int *d_neighbors,unsigned int *d_neighbors2,unsigned int *q_sig,unsigned int *d_sig,
                                        unsigned int iter,unsigned int warp_id,unsigned int *counter,unsigned int lane_id){
    unsigned int *neighbors;
    unsigned int v = pre_copy[first_index];
    unsigned int start;
    unsigned int end;
    if(join_type == 0){
        neighbors = d_neighbors;
        start = offset[v];
        end = offset[v+1];
    }else{
        neighbors = d_neighbors2;
        start = offset2[v];
        end = offset2[v+1];
    }
    for(unsigned int i=lane_id+start;i<end;i+=32){
        unsigned int temp_node = neighbors[i];
        if(compare_signature_gpu(q_sig,&d_sig[temp_node*Signature_Properties])&&!basic_search(temp_node,pre_copy,iter)){
            //if(compare_signature_gpu(q_sig,&d_sig[temp_node*Signature_Properties])&&!basic_search(temp_node,pre_copy,iter)){
            unsigned int index = atomicAdd(&counter[warp_id],1);
            if(index>=MAX_NE){
                helper_buffer[index-MAX_NE] = temp_node;
            }else{
                intersection[index] = temp_node;
            }
            intersection_bits[temp_node] = 1;
        }
    }
}
__device__ void do_intersection(unsigned int *pre_copy,unsigned int index,char join_type,char *intersection_bits,
                                unsigned int *offset,unsigned int *offset2,
                                unsigned int *d_neighbors,unsigned int *d_neighbors2,
                                unsigned int iter,unsigned int warp_id,
                                unsigned int lane_id,unsigned int t){
    unsigned int *neighbors;
    unsigned int v = pre_copy[index];
    unsigned int start;
    unsigned int end;
    if(join_type == 0){
        neighbors = d_neighbors;
        start = offset[v];
        end = offset[v+1];
    }else{
        neighbors = d_neighbors2;
        start = offset2[v];
        end = offset2[v+1];
    }
    for(unsigned int i=lane_id+start;i<end;i+=32){
        unsigned int temp_node = neighbors[i];
        if(intersection_bits[temp_node]==t&&!basic_search(temp_node,pre_copy,iter)){
            //if(intersection_bits[temp_node]==t&&!basic_search(temp_node,pre_copy,iter)){
            intersection_bits[temp_node]++;
        }
    }
}
__device__ void clean_bits(char *bits,unsigned int count,unsigned int lane_id){
    for(unsigned int i=lane_id;i<count;i+=32){
        bits[i] = 0;
    }
}
__device__ void initialize_joins(unsigned int *parents_offset,unsigned int *succs_offset,
                                 unsigned int *parents,unsigned int * succs,unsigned int v,
                                 unsigned int *joins, char * joins_type,unsigned int *join_len,unsigned int lane_id){
    unsigned int joins_count = 0;
    unsigned int start = parents_offset[v];
    unsigned int end = parents_offset[v+1];
    joins_count += (end - start);
    for(unsigned int i=lane_id;i<end - start;i+=32){
        joins[i] = parents[start + i];
        joins_type[i] = 0;
    }
    start = succs_offset[v];
    end = succs_offset[v+1];
    for(unsigned int i=lane_id;i<end - start;i+=32){
        joins[joins_count+i] = succs[start + i];
        joins_type[joins_count+i] = 1;
    }
    joins_count += (end - start);
    join_len[0] = joins_count;
}
__device__ unsigned int find_vertex_has_least_neighbor(unsigned int *offset1,unsigned int *offset2,
                                                       unsigned int *pre_copy,unsigned int *joins,
                                                       unsigned int join_len,char *joins_type){
    unsigned int least_neighbor = 10000000;
    unsigned int least_index = 0;
    for(unsigned int i=0;i<join_len;++i){
        char type = joins_type[i];
        unsigned int v = pre_copy[joins[i]];
        unsigned int neighbors_count;
        if(type == 0){
            neighbors_count = offset1[v+1] - offset1[v];
        }else{
            neighbors_count = offset2[v+1] - offset2[v];
        }
        if(least_neighbor>neighbors_count){
            least_neighbor = neighbors_count;
            least_index = i;
        }
    }
    return least_index;
}
__global__ void search_kernel(unsigned int *q_neighbors,unsigned int *q_neighbors_offset,
                              unsigned int *d_neighbors,unsigned int *d_neighbors_offset,
                              unsigned int *q_signatures,unsigned int *d_signatures,
                              unsigned int *parents,unsigned int *result_table,
                              unsigned long long int *lengths,unsigned int *parents_offset,unsigned int *succs,
                              unsigned int *succs_offset,unsigned int *order_sequence,unsigned int *reverse_neighbors,
                              unsigned int *offset2,unsigned int V,unsigned int U,
                              unsigned int iter,unsigned int *helper_buffer,char *bits_helper){
    __shared__ unsigned int signature[Signature_Properties];
    __shared__ unsigned int join_len[1];
    __shared__ unsigned int joins[2*QUERY_NODES];
    __shared__ char joins_type[2*QUERY_NODES];
    __shared__ unsigned int intersections[WARPS_EACH_BLK*MAX_NE];
    __shared__ unsigned int pre_copy[WARPS_EACH_BLK*QUERY_NODES];
    __shared__ unsigned int counter[WARPS_EACH_BLK];
    unsigned int warp_id = threadIdx.x/32;
    unsigned int lane_id = threadIdx.x%32;
    unsigned int v = order_sequence[iter];
    unsigned int intersection_len = 0;
    unsigned int global_idx = (blockIdx.x)*WARPS_EACH_BLK+warp_id;
    unsigned int helperOffset = global_idx * HelperSize;
    unsigned long long int job_from;
    unsigned long long int job_to;
    unsigned long long int job_end;
    //unsigned int mask = 0xFFFFFFFF;
    job_from = ((lengths[0]-1)/(BLK_NUMS*WARPS_EACH_BLK)+1)*global_idx;
    job_to = ((lengths[0]-1)/(BLK_NUMS*WARPS_EACH_BLK)+1)*(global_idx+1);
    job_end = lengths[0];
    if(warp_id == 0){
        initialize_joins(parents_offset,succs_offset,parents,succs,v,joins,joins_type,join_len,lane_id);
    }
    if(warp_id == 1){
        memory_copy(&q_signatures[v*Signature_Properties],signature,lane_id,Signature_Properties);
    }
    __syncthreads();
    if(job_to>=job_end){
        job_to = job_end;
    }
    for(unsigned long long int pre_idx = job_from;pre_idx<job_to;pre_idx++) {
        intersection_len = 0;
        unsigned int least_index = 0;
        for(unsigned int i = lane_id; i < iter; i += 32) {
            pre_copy[warp_id * QUERY_NODES + i] = result_table[pre_idx * (iter) + i];
        }
        if(lane_id == 0){
            counter[warp_id] = 0;
        }
        __syncwarp();
        if(lane_id == 0){
            least_index = find_vertex_has_least_neighbor(d_neighbors_offset,offset2,
                                                         &pre_copy[warp_id*QUERY_NODES],joins,join_len[0],joins_type);
        }
        least_index = __shfl_sync(0xFFFFFFFF,least_index,0);
        initialize_intersection(&pre_copy[warp_id * QUERY_NODES], joins[least_index], joins_type[least_index],
                                &intersections[warp_id * MAX_NE], &helper_buffer[helperOffset], &bits_helper[global_idx * U],
                                d_neighbors_offset, offset2,
                                d_neighbors, reverse_neighbors, signature, d_signatures, iter, warp_id,
                                counter,lane_id);
        __syncwarp();
        intersection_len = counter[warp_id];
        if(intersection_len == 0){
            continue;
        }
        __syncwarp();
        unsigned int intersection_count = 1;
        for(unsigned int i=0;i<join_len[0];i++){
            if(i == least_index){
                continue;
            }
            do_intersection(&pre_copy[warp_id * QUERY_NODES],joins[i],joins_type[i],&bits_helper[global_idx * U],
                            d_neighbors_offset,offset2,
                            d_neighbors,reverse_neighbors,
                            iter,warp_id,lane_id,intersection_count);
            intersection_count++;
            __syncwarp();
        }
        for(unsigned int i=lane_id;i<intersection_len;i+=32){
            unsigned int candidate;
            if(i<MAX_NE){
                candidate = intersections[warp_id*MAX_NE+i];
            }else{
                candidate = helper_buffer[helperOffset+i-MAX_NE];
            }
            if(bits_helper[global_idx * U+candidate] == join_len[0]){
                unsigned long long int write_offset = atomicAdd(&lengths[1],1);
                for(unsigned int j=0;j<iter;++j){
                    result_table[lengths[0]*iter+(write_offset)*(iter+1)+j] = pre_copy[warp_id * QUERY_NODES + j];
                }
                result_table[lengths[0]*iter+(write_offset)*(iter+1)+iter] = candidate;
            }
        }
        __syncwarp();
    }
}
unsigned int search(Graph &query_graph, Graph &data_graph, string queryFile,string dataFile,unsigned int * result_table){
    string queryBaseName = getFileBaseName(queryFile);
    string dataBaseName = getFileBaseName(dataFile);
    unsigned int iters = query_graph.V;
    unsigned int U = data_graph.V;
    unsigned int *order_sequence = query_graph.order_sequence;
    unsigned int *parents = query_graph.parents;
    unsigned int *parents_offset = query_graph.parents_offset;
    unsigned int *succs = query_graph.succs;
    unsigned int *succs_offset = query_graph.succs_offset;
    cudaEvent_t event_start;
    cudaEvent_t event_stop;
    cudaEventCreate(&event_start);
    cudaEventCreate(&event_stop);
    //declare gpu memory
    unsigned int *order_sequence_gpu;
    unsigned int *parents_gpu;
    unsigned int *parents_gpu_offset;
    unsigned int *succs_gpu;
    unsigned int *succs_gpu_offset;
    unsigned int *results;
    unsigned long long int *lengths;
    unsigned int *q_neighbors_gpu;
    unsigned int *q_neighbors_gpu_offset;
    unsigned int *d_neighbors_gpu;
    unsigned int *d_neighbors_gpu_offset;
    unsigned int *q_signatures;
    unsigned int *d_signatures;
    unsigned int *reverse_neighbors;
    unsigned int *offset2;
    unsigned int *helperBuffer;
    char * bits_helper;
    dim3 grid(BLK_NUMS);
    dim3 block(BLK_DIM);
    //malloc gpu memory
    cudaEventRecord(event_start);
    cudaMalloc(&order_sequence_gpu,iters*sizeof(unsigned int));
    cudaMalloc(&parents_gpu,parents_offset[iters]*sizeof(unsigned int));
    cudaMalloc(&parents_gpu_offset,(iters+1)*sizeof(unsigned int));
    cudaMalloc(&succs_gpu,succs_offset[iters]*sizeof(unsigned int));
    cudaMalloc(&succs_gpu_offset,(iters+1)*sizeof(unsigned int));
    cudaMalloc(&q_neighbors_gpu,query_graph.offset[iters]*sizeof(unsigned int));
    cudaMalloc(&q_neighbors_gpu_offset,(iters+1)*sizeof(unsigned int));
    cudaMalloc(&d_neighbors_gpu,data_graph.offset[U]*sizeof(unsigned int));
    cudaMalloc(&d_neighbors_gpu_offset,(U+1)*sizeof(unsigned int));
    cudaMalloc(&reverse_neighbors,data_graph.offset2[U]*sizeof(unsigned int));
    cudaMalloc(&offset2,(U+1)*sizeof(unsigned int));
    cudaMalloc(&q_signatures,iters*Signature_Properties*sizeof(unsigned int));
    cudaMalloc(&d_signatures,U*Signature_Properties*sizeof(unsigned int));
    cudaMalloc(&bits_helper,WARPS_EACH_BLK*BLK_NUMS*U*sizeof(char));
    cudaMalloc(&results,GPU_TABLE_LIMIT*sizeof(unsigned int));
    cudaMalloc(&lengths,2*sizeof(unsigned long long int));
    cudaMalloc(&helperBuffer,WARPS_EACH_BLK*BLK_NUMS*HelperSize*sizeof(unsigned int));
    chkerr(cudaDeviceSynchronize());
    cudaEventRecord(event_stop);
    cudaEventSynchronize(event_stop);
    float time_cuda_malloc = 0;
    cudaEventElapsedTime(&time_cuda_malloc, event_start, event_stop);
    cudaEventRecord(event_start);
    cudaMemcpy(order_sequence_gpu,order_sequence,iters*sizeof(unsigned int),cudaMemcpyHostToDevice);
    cudaMemcpy(parents_gpu,parents,parents_offset[iters]*sizeof(unsigned int),cudaMemcpyHostToDevice);
    cudaMemcpy(parents_gpu_offset,parents_offset,(iters+1)*sizeof(unsigned int),cudaMemcpyHostToDevice);
    cudaMemcpy(succs_gpu,succs,succs_offset[iters]*sizeof(unsigned int),cudaMemcpyHostToDevice);
    cudaMemcpy(succs_gpu_offset,succs_offset,(iters+1)*sizeof(unsigned int),cudaMemcpyHostToDevice);
    cudaMemset(lengths,0,(2)*sizeof(unsigned long long int));
    cudaMemcpy(q_neighbors_gpu,query_graph.neighbors,query_graph.offset[iters]*sizeof(unsigned int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(q_neighbors_gpu_offset,query_graph.offset,(iters+1)*sizeof(unsigned int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighbors_gpu,data_graph.neighbors,data_graph.offset[U]*sizeof(unsigned int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_neighbors_gpu_offset,data_graph.offset,(U+1)*sizeof(unsigned int),cudaMemcpyHostToDevice);
    cudaMemcpy(q_signatures,query_graph.signatures,iters*Signature_Properties*sizeof(unsigned int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_signatures,data_graph.signatures,U*Signature_Properties*sizeof(unsigned int),cudaMemcpyHostToDevice);
    cudaMemcpy(q_signatures,query_graph.signatures,iters*Signature_Properties*sizeof(unsigned int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(reverse_neighbors,data_graph.reverse_neighbors,data_graph.offset2[U]*sizeof(unsigned int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(offset2,data_graph.offset2,(U+1)*sizeof(unsigned int),cudaMemcpyHostToDevice);
    cudaMemset(bits_helper,0,WARPS_EACH_BLK*BLK_NUMS*U*sizeof(char));
    cudaEventRecord(event_stop);
    cudaEventSynchronize(event_stop);
    unsigned long long int results_count = 1;
    unsigned long long int startEnd[2];
    unsigned int ini_count = 0;
    cudaEventRecord(event_start);
    initialize_searching<<<68,1024>>>(q_signatures,d_signatures,results,order_sequence_gpu,data_graph.V,lengths);
    cudaMemcpy(&ini_count,&lengths[0],1*sizeof(unsigned long long int),cudaMemcpyDeviceToHost);
    unsigned int *cpu_buffer = (unsigned int *)malloc(U*sizeof(unsigned int));
    chkerr(cudaMemcpy(cpu_buffer,results,ini_count*sizeof(unsigned int),cudaMemcpyDeviceToHost));
    shuffle_array(cpu_buffer,ini_count);
    unsigned int final_count = 0;
    unsigned int t_size = 512;
    unsigned int loops = (ini_count -1)/t_size + 1;
    //cout<<ini_count<<endl;
    if(ini_count>0){
        for(unsigned int l=0;l<loops;l++){
            cudaMemset(bits_helper,0,WARPS_EACH_BLK*BLK_NUMS*U*sizeof(char));
            startEnd[0] = t_size;
            if(startEnd[0]>ini_count - l*t_size){
                startEnd[0] = ini_count - l*t_size;
            }
            chkerr(cudaMemset(lengths,0,2*sizeof(unsigned int)));
            chkerr(cudaMemcpy(lengths,startEnd,1*sizeof(unsigned int),cudaMemcpyHostToDevice));
            chkerr(cudaMemcpy(results,&cpu_buffer[l*t_size],startEnd[0]*sizeof(unsigned int),cudaMemcpyHostToDevice));
            for(unsigned int iter=1;iter<iters;++iter){
                search_kernel<<<grid,block>>>(q_neighbors_gpu,q_neighbors_gpu_offset,
                        d_neighbors_gpu,d_neighbors_gpu_offset,
                        q_signatures,d_signatures,
                        parents_gpu,results,
                        lengths,parents_gpu_offset,succs_gpu,
                        succs_gpu_offset,order_sequence_gpu,
                        reverse_neighbors,offset2,
                        query_graph.V,data_graph.V,
                        iter,helperBuffer,bits_helper);
                chkerr(cudaDeviceSynchronize());
                cudaMemcpy(startEnd,&lengths[0],2*sizeof(unsigned long long int),cudaMemcpyDeviceToHost);
                results_count = startEnd[1];
                if(results_count == 0){
                    break;
                }
                unsigned int preBufferSize = startEnd[0]*iter;
                unsigned int currBufferSize = startEnd[1]*(iter+1);
                unsigned int copy_iters = (currBufferSize - 1)/preBufferSize + 1;
                for(unsigned int copy_iter = 0;copy_iter<copy_iters;++copy_iter){
                    unsigned int trunk_size = preBufferSize;
                    if(trunk_size > (currBufferSize - copy_iter*preBufferSize)){
                        trunk_size = currBufferSize - copy_iter*preBufferSize;
                    }
                    cudaMemcpy(&results[copy_iter*preBufferSize],&results[preBufferSize+copy_iter*preBufferSize],trunk_size*sizeof(unsigned int),
                               cudaMemcpyDeviceToDevice);
                }
                cudaMemcpy(&lengths[0],&lengths[1],1*sizeof(unsigned long long int),cudaMemcpyDeviceToDevice);
                cudaMemset(&lengths[1],0,1*sizeof(unsigned long long int));
            }
            if(results_count>0){
                cudaMemcpy(result_table,results,results_count*query_graph.V,cudaMemcpyDeviceToHost);
                final_count += results_count;
            }
        }
    }
    cudaEventRecord(event_stop);
    cudaEventSynchronize(event_stop);
    float time_milli_sec = 0;
    cudaEventElapsedTime(&time_milli_sec, event_start, event_stop);
    cout<<queryBaseName<<","<<dataBaseName<<","<<final_count<<","<<time_milli_sec<<endl;
    return final_count;
}
int main(int argc, char *argv[]){
    std::string query_graph_file = argv[2];
    std::string data_graph_file = argv[1];
    Graph query_graph(0,query_graph_file);
    Graph data_graph(1,data_graph_file);
    unsigned int *result_table = (unsigned int *)malloc(GPU_TABLE_LIMIT*sizeof(unsigned int));
    unsigned int result_len = search(query_graph,data_graph,query_graph_file,data_graph_file,result_table);
    return 0;
}

