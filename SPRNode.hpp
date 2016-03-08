#ifndef SPRNode_hpp
#define SPRNode_hpp

#define COPY_CONTRACTED

#include <stdio.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <set>
#include "Forest.hpp"
struct StringCompare {
    bool operator() (const std::string &a, const std::string &b) const {
        return strcmp(a.c_str(), b.c_str()) < 0;
    }
};


// representation of a component with no leaves
#define DEAD_COMPONENT "*"
//#define DEBUG_PROTECTED "@"
/*void find_sibling_pairs_hlpr(SPRNode *SPRNode, list<SPRNode *> &sibling_pairs);
 void find_leaves_hlpr(SPRNode *SPRNode, vector<SPRNode *> &leaves);
 void str_subtree_hlpr(string *s);
 list<SPRNode *> find_sibling_pairs(SPRNode *SPRNode);
 vector<SPRNode *> find_leaves(SPRNode *SPRNode);
 */
int stomini(std::string s);

class Forest;

class SPRNode {
private:
    //SPRNode *lc;			// left child
    //SPRNode *rc;			// right child
    std::list<SPRNode *> children; // children
    SPRNode *p;			// parent
    std::list<SPRNode *>:: iterator p_link;		// location in parents list
    SPRNode *twin;			// counterpart in another tree
    std::string name;		// label
    int SPR_depth;			//distance from root
    int pre_num;	// preorder number
    int edge_pre_start;
    int edge_pre_end;
    
    int component_number;
    std::list <SPRNode *> active_descendants;
    std::list <SPRNode *> root_lcas;
    std::list<std::list <SPRNode *>::iterator> removable_descendants;
    std::list <SPRNode *>::iterator sibling_pair_loc;
    int sibling_pair_status;
    int num_clustered_children;
    Forest *forest;
    // TODO: contracted_list ?
    SPRNode *contracted_lc;
    SPRNode *contracted_rc;
    bool contracted;
    bool edge_protected;
    int max_merge_SPR_depth;
    bool allow_sibling;
    int lost_children;
    double support;
    double support_normalization;

public:
    SPRNode();
    SPRNode(std::string n);
    SPRNode(std::string n, int d);
    SPRNode(SPRNode *lc, SPRNode *rc, SPRNode *p, std::string n, int d);
    void init(SPRNode *lc, SPRNode *rc, SPRNode *p, std::string n, int d);
    SPRNode(const SPRNode &n);
    SPRNode(const SPRNode &n, SPRNode *parent);
    ~SPRNode();
    void delete_tree();
    void delete_child(SPRNode *n);
    void add_child(SPRNode *n);
    void add_child_keep_SPR_depth(SPRNode *n);
    void insert_child(SPRNode *sibling, SPRNode *n);
    void insert_child_keep_SPR_depth(SPRNode *sibling, SPRNode *n);
    SPRNode *set_twin(SPRNode *n);
    void set_name(std::string n);
    int set_SPR_depth(int d);
    void fix_SPR_depths();
    int set_preorder_number(int p);
    int set_edge_pre_start(int p);
    int set_edge_pre_end(int p);
    void copy_edge_pre_interval(SPRNode *n);
    int set_component_number(int c);
    std::list<SPRNode *>& get_children();
    SPRNode *get_contracted_lc();
    SPRNode *get_contracted_rc();
    SPRNode *set_contracted_lc(SPRNode *n);
    SPRNode *set_contracted_rc(SPRNode *n);
    bool is_protected();
    bool is_contracted();
    void protect_edge();
    void unprotect_edge();
    void unprotect_subtree();
    void protect_supported_edges();
    bool can_be_sibling();
    void disallow_siblings();
    void disallow_siblings_subtree();
    void allow_siblings();
    void allow_siblings_subtree();
    int num_lost_children();
    int get_max_merge_SPR_depth();
    void set_max_merge_SPR_depth(int d);
    int count_lost_children_subtree();
    int count_lost_subtree();
    void lost_child();
    void no_lost_children();
    void no_lost_children_subtree();
    double get_support();
    void set_support(double s);
    void a_inc_support();
    void a_dec_support();
    double get_support_normalization();
    void set_support_normalization(double s);
    void a_inc_support_normalization();
    void a_dec_support_normalization();
    void normalize_support();
    int get_component_number();
    void increase_clustered_children();
    void decrease_clustered_children();
    void set_num_clustered_children(int c);
    int get_num_clustered_children();
    std::list <SPRNode *> *get_active_descendants();
    std::list <SPRNode *> *get_root_lcas();
    int get_sibling_pair_status();
    void set_sibling_pair_status(int s);
    void set_forest(Forest *f);
    void set_forest_rec(Forest *f);
    Forest *get_forest();
    std::list<std::list <SPRNode *>::iterator> *get_removable_descendants();
    void initialize_component_number(int value);
    void initialize_active_descendants(std::list <SPRNode *> value);
    void initialize_root_lcas(std::list <SPRNode *> value);
    void initialize_removable_descendants(std::list<std::list <SPRNode *>::iterator> value);
    SPRNode *contract(bool remove);
    SPRNode *contract();
    bool contract_sibling_pair();
    bool contract_sibling_pair_undoable();
    SPRNode *contract_sibling_pair_undoable(SPRNode *child1, SPRNode *child2);
    void undo_contract_sibling_pair();
    void fix_contracted_order();
    void cut_parent();
    void contract_SPRNode();
    SPRNode *parent();
    inline SPRNode *lchild();
    SPRNode *rchild();
    SPRNode *get_twin();
    int get_SPR_depth();
    int get_preorder_number();
    int get_edge_pre_start();
    int get_edge_pre_end();
    std::string str();
    std::string get_name();
    void str_hlpr(std::string *s);
    std::string str_subtree();
    void str_subtree_hlpr(std::string *s);
    std::string str_support_subtree(bool allow_negative);
    std::string str_support_subtree();
    void str_support_subtree_hlpr(std::string *s, bool allow_negative);
    std::string str_edge_pre_interval_subtree();
    void str_edge_pre_interval_subtree_hlpr(std::string *s);
    void str_c_subtree_hlpr(std::string *s);
    std::string str_subtree_twin();
    void str_subtree_twin_hlpr(std::string *s);
    void print();
    void print_subtree();
    void print_subtree_hlpr();
    void print_subtree_twin_hlpr();
    bool is_leaf();
    bool is_sibling_pair();
    bool is_singleton();
    void find_sibling_pairs_hlpr(std::list<SPRNode *> *sibling_pairs);
    void append_sibling_pairs(std::list<SPRNode *> *sibling_pairs);
    std::list<SPRNode *> *find_sibling_pairs();
    void find_leaves_hlpr(std::vector<SPRNode *> &leaves);
    std::vector<SPRNode *> find_leaves();
    void find_interior_hlpr(std::vector<SPRNode *> &interior);
    std::vector<SPRNode *> find_interior();
    void find_descendants_hlpr(std::vector<SPRNode *> &descendants);
    std::vector<SPRNode *> find_descendants();
    bool contains_leaf(int number);
    void resync();
    void unsync();
    void unsync_interior();
    void sync_af_twins();
    SPRNode *find_root();
    bool same_component(SPRNode *n, int &lca_SPR_depth, int &path_length);
    bool same_component(SPRNode *n);
    bool same_component(SPRNode *n, int &lca_SPR_depth);
    void labels_to_numbers(std::map<std::string, int> *label_map, std::map<int, std::string> *reverse_label_map);
    void numbers_to_labels(std::map<int, std::string> *reverse_label_map);
    void build_name_to_pre_map(std::map<std::string, int> *name_to_pre);
    void count_numbered_labels(std::vector<int> *label_counts);
    void preorder_number();
    int preorder_number(int next);
    void edge_preorder_interval();
    int size();
    int size_using_prenum();
    int max_SPR_depth();
    size_t max_degree();
    void add_to_front_sibling_pairs(std::list<SPRNode *> *sibling_pairs, int status);
    void add_to_sibling_pairs(std::list<SPRNode *> *sibling_pairs, int status);
    void remove_sibling_pair(std::list<SPRNode *> *sibling_pairs);
    void clear_sibling_pair(std::list<SPRNode *> *sibling_pairs);
    SPRNode *get_sibling();
    SPRNode *get_right_sibling();
    SPRNode *get_sibling(std::list<SPRNode *> *sibling_pairs);
    void set_sibling(SPRNode *sibling);
    void clear_sibling_pair_status();
    void fix_parents();
    void left_rotate();
    void right_rotate();
    void next_rooting();
    void reroot(SPRNode *new_lc);
    void fixroot();
    SPRNode *expand_parent_edge(SPRNode *n);
    SPRNode *undo_expand_parent_edge();
    SPRNode *spr(SPRNode *new_sibling, int &which_child);
    void spr(SPRNode *new_sibling);
    void find_descendant_counts_hlpr(std::vector<int> *dc);
    std::vector<int> *find_descendant_counts();
    void find_leaf_counts_hlpr(std::vector<int> *lc);
    std::vector<int> *find_leaf_counts();
    SPRNode *find_median();
    SPRNode *find_subtree_of_size(double percentage);
    SPRNode *find_median_hlpr(std::vector<int> *dc, int target_size);
    int any_leaf_preorder_number();
    SPRNode *find_by_prenum(int prenum);
    void expand_contracted_SPRNodes();
    int get_name_num();
};

// function prototypes
SPRNode *spr_building_tree(std::string s);
SPRNode *spr_building_tree(std::string s, std::set<std::string, StringCompare> *include_only);
SPRNode *spr_building_tree(std::string s, int start_SPR_depth);
SPRNode *spr_building_tree(std::string s, int start_SPR_depth, std::set<std::string, StringCompare> *include_only);
size_t spr_building_tree_helper(size_t start, const std::string& s, SPRNode *parent,
                                bool &valid, std::set<std::string, StringCompare> *include_only);
void swap(SPRNode **a, SPRNode **b);
int stomini(std::string s);
std::string SPR_root(std::string s);

















#endif /* SPRNode_hpp */
