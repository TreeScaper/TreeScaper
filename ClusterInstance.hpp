

#ifndef ClusterInstance_hpp
#define ClusterInstance_hpp

#include "Forest.hpp"
#include "SPRNode.hpp"

class ClusterInstance{
public:
    Forest *F1;
    Forest *F2;
    SPRNode *F1_cluster_SPRNode;
    SPRNode *F2_cluster_SPRNode;
    bool F2_has_component_zero;
    
public:
    ClusterInstance();
    ClusterInstance(Forest *f1, Forest *f2, SPRNode *f1_cluster_SPRNode,
                    SPRNode *f2_cluster_SPRNode,bool f2_has_component_zero);
    ClusterInstance(const ClusterInstance &c);
    ~ClusterInstance();
    void init(Forest *f1, Forest *f2, SPRNode *f1_cluster_SPRNode,
              SPRNode *f2_cluster_SPRNode, bool f2_has_component_zero);
    bool is_original();
    int join_cluster(Forest *original_F1, Forest *original_F2);

};

void cluster_reduction_find_components(SPRNode *n,
                                       std::vector<bool> *F2_cluster_copy_components,
                                       std::vector<bool> *old_F2_keep_components,
                                       int cluster_component_number);
std::list<ClusterInstance> cluster_reduction(Forest *old_F1, Forest *old_F2,
                                             std::list<SPRNode *> *cluster_points);

#endif /* ClusterInstance_hpp */
