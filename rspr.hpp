#ifndef rspr_hpp
#define rspr_hpp

#include <stdio.h>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <set>
#include <list>
#include "Forest.hpp"
#include "ClusterForest.hpp"
#include "SPRLCA.hpp"
#include "ClusterInstance.hpp"
#include "SiblingPair.hpp"
#include "UndoMachine.hpp"

using namespace std;


enum RELAXATION {SPR_STRICT, NEGATIVE_RELAXED, ALL_RELAXED};

int rSPR_3_approx_hlpr(Forest *T1, Forest *T2, std::list<SPRNode *> *singletons,
                       std::list<SPRNode *> *sibling_pairs);
int rSPR_3_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_approx_hlpr(Forest *T1, Forest *T2, std::list<SPRNode *> *singletons, std::list<SPRNode *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests);
int rSPR_worse_3_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_approx(Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx(SPRNode *subtree, Forest *T1, Forest *T2);
int rSPR_worse_3_approx(SPRNode *subtree, Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx_binary_hlpr(Forest *T1, Forest *T2, std::list<SPRNode *> *singletons, std::list<SPRNode *> *sibling_pairs, Forest **F1, Forest **F2, bool save_forests);
int rSPR_worse_3_approx_binary(Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx_binary(Forest *T1, Forest *T2);
int rSPR_branch_and_bound(Forest *T1, Forest *T2);
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k);
int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int end_k);
int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int start_k, int end_k);
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
                               std::set<SiblingPair> *sibling_pairs, std::list<SPRNode *> *singletons, bool cut_b_only,
                               std::list<std::pair<Forest,Forest> > *AFs, std::list<SPRNode *> *protected_stack,
                               int *num_ties);
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
                               std::set<SiblingPair> *sibling_pairs, std::list<SPRNode *> *singletons, bool cut_b_only,
                               std::list<std::pair<Forest,Forest> > *AFs, std::list<SPRNode *> *protected_stack,
                               int *num_ties, SPRNode *prev_T1_a, SPRNode *prev_T1_c);
int rSPR_total_approx_distance(SPRNode *T1, std::vector<SPRNode *> &gene_trees);
int rSPR_total_approx_distance(SPRNode *T1, std::vector<SPRNode *> &gene_trees,
                               int threshold);
int rSPR_total_distance(SPRNode *T1, std::vector<SPRNode *> &gene_trees);
int rSPR_total_distance(SPRNode *T1, std::vector<SPRNode *> &gene_trees,
                        std::vector<int> *original_scores);

void rSPR_pairwise_distance(SPRNode *T1, std::vector<SPRNode *> &gene_trees);
void rSPR_pairwise_distance(SPRNode *T1, std::vector<SPRNode *> &gene_trees, bool approx);
void rSPR_pairwise_distance(SPRNode *T1, std::vector<SPRNode *> &gene_trees, int start, int end);
void rSPR_pairwise_distance(SPRNode *T1, std::vector<SPRNode *> &gene_trees, int start, int end, bool approx);
void rSPR_pairwise_distance(SPRNode *T1, std::vector<SPRNode *> &gene_trees, int max_spr);
void rSPR_pairwise_distance(SPRNode *T1, std::vector<SPRNode *> &gene_trees, int max_spr, int start, int end);
void rSPR_pairwise_distance_unrooted(SPRNode *T1, std::vector<SPRNode *> &gene_trees);
void rSPR_pairwise_distance_unrooted(SPRNode *T1, std::vector<SPRNode *> &gene_trees, int start, int end);
void rSPR_pairwise_distance_unrooted(SPRNode *T1, std::vector<SPRNode *> &gene_trees, bool approx);
void rSPR_pairwise_distance_unrooted(SPRNode *T1, std::vector<SPRNode *> &gene_trees, int start, int end, bool approx);

int rSPR_total_distance_unrooted(SPRNode *T1, std::vector<SPRNode *> &gene_trees, int threshold);
int rSPR_total_distance_unrooted(SPRNode *T1, std::vector<SPRNode *> &gene_trees, int threshold, std::vector<int> *original_scores);

int rSPR_branch_and_bound_simple_clustering(SPRNode *T1, SPRNode *T2, Forest **out_F1, Forest **out_F2);
int rSPR_branch_and_bound_simple_clustering(SPRNode *T1, SPRNode *T2, bool verbose, std::map<std::string, int> *label_map, std::map<int, std::string> *reverse_label_map);
int rSPR_branch_and_bound_simple_clustering(SPRNode *T1, SPRNode *T2, bool verbose, std::map<std::string, int> *label_map, std::map<int, std::string> *reverse_label_map, int min_k, int max_k);
int rSPR_branch_and_bound_simple_clustering(SPRNode *T1, SPRNode *T2);
int rSPR_branch_and_bound_simple_clustering(SPRNode *T1, SPRNode *T2, bool verbose);
int rSPR_branch_and_bound_simple_clustering(SPRNode *T1, SPRNode *T2, bool verbose, int min_k, int max_k);
int rSPR_branch_and_bound_simple_clustering(SPRNode *T1, SPRNode *T2, bool verbose, std::map<std::string, int> *label_map, std::map<int, std::string> *reverse_label_map, int min_k, int max_k, Forest **out_F1, Forest **out_F2);

void reduction_leaf(Forest *T1, Forest *T2);
void reduction_leaf(Forest *T1, Forest *T2, UndoMachine *um);
bool chain_match(SPRNode *T1_SPRNode, SPRNode *T2_SPRNode, SPRNode *T2_SPRNode_end);
bool is_nonbranching(Forest *T1, Forest *T2, SPRNode *T1_a, SPRNode *T1_c, SPRNode *T2_a, SPRNode *T2_c);

SPRNode *find_subtree_of_approx_distance(SPRNode *n, Forest *F1, Forest *F2, int target_size);
SPRNode *find_subtree_of_approx_distance_hlpr(SPRNode *n, Forest *F1, Forest *F2, int target_size);
SPRNode *find_best_root(SPRNode *T1, SPRNode *T2);
double find_best_root_acc(SPRNode *T1, SPRNode *T2);
void find_best_root_hlpr(SPRNode *T2, int pre_separator, int group_1_total,
                         int group_2_total, SPRNode **best_root, double *best_root_b_acc);
void find_best_root_hlpr(SPRNode *n, int pre_separator, int group_1_total,
                         int group_2_total, SPRNode **best_root, double *best_root_b_acc,
                         int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties);
bool contains_bipartition(SPRNode *n, int pre_start, int pre_end,
                          int group_1_total, int group_2_total, int *p_group_1_descendants,
                          int *p_group_2_descendants);
bool outgroup_root(SPRNode *T, std::set<std::string, StringCompare> outgroup);
bool outgroup_root(SPRNode *n, std::vector<int> &num_in, std::vector<int> &num_out);
bool outgroup_reroot(SPRNode *n, std::vector<int> &num_in, std::vector<int> &num_out);
void count_in_out(SPRNode *n, std::vector<int> &num_in, std::vector<int> &num_out,
                  std::set<std::string, StringCompare> &outgroup);


/*Joel's part*/
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose, std::map<std::string, int> *label_map,std::map<int, std::string> *reverse_label_map);
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2);
int rSPR_branch_and_bound_simple_clustering(Forest *T1, Forest *T2, bool verbose);
int rSPR_total_distance(Forest *T1, std::vector<SPRNode *> &gene_trees);

int rSPR_worse_3_approx_distance_only(Forest *T1, Forest *T2);


class ProblemSolution {
public:
    std::string T1;
    std::string T2;
    int k;
    
    ProblemSolution(Forest *t1, Forest *t2, int new_k)
    {
        T1 = t1->str();
        T2 = t2->str();
        k = new_k;
    }
};

void add_sibling_pair(std::set<SiblingPair> *sibling_pairs, SPRNode *a, SPRNode *c, UndoMachine *um);
SiblingPair pop_sibling_pair(std::set<SiblingPair> *sibling_pairs, UndoMachine *um);
SiblingPair pop_sibling_pair(std::set<SiblingPair>::iterator s, std::set<SiblingPair> *sibling_pairs, UndoMachine *um);
bool chain_match(SPRNode *T1_SPRNode, SPRNode *T2_SPRNode, SPRNode *T2_SPRNode_end);



























#endif /* rspr_hpp */
