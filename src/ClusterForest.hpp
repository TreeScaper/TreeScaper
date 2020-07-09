#ifndef ClusterForest_hpp
#define ClusterForest_hpp

#include "Forest.hpp"

class ClusterForest: public Forest{
public:
    std::vector<SPRNode *> cluster_nodes;

public:
    ClusterForest();
    ClusterForest(std::vector<SPRNode *> components);
    ClusterForest(SPRNode *head);
    ClusterForest(const ClusterForest &f);
    void init();
    ~ClusterForest();
    
    void add_cluster(SPRNode *n, std::string name);
    // ClusterForest_ClusterForest_swap the contents of two forests
    void ClusterForest_swap(ClusterForest *f);
    void ClusterForest_swap(Forest *f);
    SPRNode *get_cluster_node(int i);
    size_t spr_num_clusters();
    void join_cluster(int cluster_loc, Forest *solved_cluster);
    void join_cluster(Forest *solved_cluster);
};
    
    // ClusterForest_swap two cluster forests
    void ClusterForest_swap(ClusterForest **a, ClusterForest **b);

#endif /* ClusterForest_hpp */
