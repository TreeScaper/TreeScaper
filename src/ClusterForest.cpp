#include "ClusterForest.hpp"

ClusterForest::ClusterForest():Forest()
{
    init();
}

ClusterForest::ClusterForest(std::vector<SPRNode *> components):Forest(components)
{
    init();
}

ClusterForest::ClusterForest(SPRNode *head):Forest(head)
{
    init();
}

ClusterForest::ClusterForest(const ClusterForest &f):Forest(f)
{
    cluster_nodes = std::vector<SPRNode *>(f.cluster_nodes);
}

void ClusterForest::init()
{
    cluster_nodes = std::vector<SPRNode *>();
    cluster_nodes.push_back(NULL);
}

ClusterForest::~ClusterForest()
{
    /* don't delete cluster SPRNodes because they are, by
		   definition, in the components
     */
    cluster_nodes.clear();
}

void ClusterForest::add_cluster(SPRNode *n, std::string name) {
    SPRNode *n_parent = n->parent();
    n->cut_parent();
    SPRNode *n_child = new SPRNode(name);
    n_parent->add_child(n_child);
    add_component(n);
    cluster_nodes.push_back(n_child);
}

// ClusterForest_ClusterForest_swap the contents of two forests
void ClusterForest::ClusterForest_swap(ClusterForest *f) {
    std::vector<SPRNode *> components_temp = this->components;
    this->components = f->components;
    f->components = components_temp;
    
    /*
     vector<SPRNode *> deleted_SPRNodes_temp = this->deleted_SPRNodes;
     this->deleted_SPRNodes = f->deleted_SPRNodes;
     f->deleted_SPRNodes = deleted_SPRNodes_temp;
     */
    
    std::vector<SPRNode *> cluster_nodes_temp = this->cluster_nodes;
    this->cluster_nodes = f->cluster_nodes;
    f->cluster_nodes = cluster_nodes_temp;
}

void ClusterForest::ClusterForest_swap(Forest *f) {
    std::vector<SPRNode *> components_temp = this->components;
    this->components = f->components;
    f->components = components_temp;
    
    if (contains_rho())
        f->rho = true;
}

SPRNode *ClusterForest::get_cluster_node(int i) {
    return cluster_nodes[i];
}

size_t ClusterForest::spr_num_clusters() {
    return cluster_nodes.size();
}

void ClusterForest::join_cluster(int cluster_loc, Forest *solved_cluster) {
    SPRNode *cluster_node = get_cluster_node(cluster_loc);
    SPRNode *cluster_parent = cluster_node->parent();
    cluster_node->cut_parent();
    //		delete cluster_node;
    cluster_node->delete_tree();
    //		delete components[cluster_loc];
    components[cluster_loc]->delete_tree();
    components[cluster_loc] = NULL;
    int start = 0;
    if (solved_cluster->contains_rho()) {
        if (cluster_parent->parent() != NULL ||
            !((cluster_parent->lchild() && cluster_parent->lchild()->get_name() == "X")
              || (cluster_parent->rchild() && cluster_parent->rchild()->get_name() == "X")))
            cluster_parent->contract(true);
    }
    else {
        cluster_parent->add_child(solved_cluster->get_component(0));
        //cluster_parent->add_child(new SPRNode(*(solved_cluster->get_component(0))));
        start = 1;
        if (cluster_parent->get_children().size() == 1) {
            cluster_parent->contract(true);
        }
    }
    // should we add these to a finished_components or something?
    for(int i = start; i < solved_cluster->num_components(); i++) {
        if (solved_cluster->get_component(i)->str() != "p")
            add_component(solved_cluster->get_component(i));
        //add_component(new SPRNode(*(solved_cluster->get_component(i))));
        else
            solved_cluster->get_component(i)->delete_tree();
    }
    solved_cluster->erase_components();
    
}

void ClusterForest::join_cluster(Forest *solved_cluster) {
    int start = 0;
    if (solved_cluster->contains_rho()) {
        add_rho();
    }
    else {
        add_component(0,solved_cluster->get_component(0));
        //add_component(0,new SPRNode(*(solved_cluster->get_component(0))));
        start = 1;
    }
    // should we add these to a finished_components or something?
    for(int i = start; i < solved_cluster->num_components(); i++) {
        if (solved_cluster->get_component(i)->str() != "p")
            //add_component(new SPRNode(*(solved_cluster->get_component(i))));
            add_component(solved_cluster->get_component(i));
        else
            solved_cluster->get_component(i)->delete_tree();
    }
    solved_cluster->erase_components();
}

// ClusterForest_swap two cluster forests
void ClusterForest_swap(ClusterForest **a, ClusterForest **b) {
    (*a)->ClusterForest_swap(*b);
}


