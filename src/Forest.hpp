#ifndef Forest_hpp
#define Forest_hpp

#include <stdio.h>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <list>
#include <deque>
#include "SPRLCA.hpp"
#include "SPRNode.hpp"
#include <map>
#include <limits>
#include <climits>


class ClusterInstance;
class SPRNode;
class LCA;

class Forest{
public:
    std::vector<SPRNode *> components;
    std::vector<SPRNode *> deleted_nodes;
    bool rho;
    Forest *twin;
    ClusterInstance *cluster;

public:
    Forest();
    Forest(std::vector<SPRNode *> components);
    Forest(SPRNode *head);
    Forest(const Forest &f);
    Forest(Forest *f);
    Forest(Forest *f, bool b);
    void init(std::vector<SPRNode *> components);
    ~Forest();
    void swap(Forest *f);
    void print_components();
    void print_components(std::string s);
    void print_components_with_twins();
    void print_components_with_edge_pre_interval();
    std::string str();
    void add_component(SPRNode *head);
    void add_component(int pos, SPRNode *head);
    void add_deleted_node(SPRNode *n);
    SPRNode *get_component(int i);
    void set_component(int i, SPRNode *head);
    inline size_t num_components();
    void set_twin(Forest *f);
    Forest *get_twin();
    ClusterInstance *get_cluster();
    void set_cluster(ClusterInstance *c);
    void set_cluster(ClusterInstance &c);
    void update_component(SPRNode *old_c, SPRNode *new_c);
    std::list<SPRNode *> *find_sibling_pairs();
    std::list<SPRNode *> find_singletons();
    void resync();
    void unsync();
    void unsync_interior();
    SPRNode *find_by_prenum(int prenum);
    void labels_to_numbers(std::map<std::string, int> *label_map, std::map<int, std::string> *reverse_label_map);
    void numbers_to_labels(std::map<int, std::string> *reverse_label_map);
    size_t size();
    bool add_rho();
    void set_rho(bool b);
    bool contains_rho();
    void erase_components(int start, int end);
    void erase_components();
    void label_SPRNodes_with_forest();
    void move_first_component_to_end() ;
    void unprotect_edges();

};

std::vector<SPRNode *> find_labels(std::vector<SPRNode *> components);
bool sync_twins(Forest *T1, Forest *T2);
void sync_interior_twins(Forest *T1, Forest *T2);
void sync_interior_twins(SPRNode *n, LCA *twin_LCA);
void sync_interior_twins(SPRNode *n, std::vector<LCA> *F2_LCAs);
std::list<SPRNode *> *find_cluster_points(Forest *F1, Forest *F2);
void find_cluster_points(SPRNode *n, std::list<SPRNode *> *cluster_points,
                         std::vector<int> *leaf_counts_F1, std::vector<int> *leaf_counts_F2);
void delete_and_merge_LCAs(std::list<SPRNode *> *active_descendants,
                           std::vector<LCA> *F2_LCAs, std::list<SPRNode *>:: iterator SPRNode1_location,
                           std::list<SPRNode *>:: iterator SPRNode2_location);
void delete_and_merge_LCAs(SPRNode *n, std::list<SPRNode *> *active_descendants,
                           std::vector<LCA> *F2_LCAs);


void sync_interior_twins_real(Forest *T1, Forest *F2);
void sync_af_twins(Forest *F1, Forest *F2);
void swap(Forest **a, Forest **b);
void expand_contracted_SPRNodes(Forest *F);
Forest *build_finished_forest(std::string &name);
Forest *build_forest(std::string &name);


#endif /* Forest_hpp */
