#include "ClusterInstance.hpp"

ClusterInstance::ClusterInstance()
{
    init(NULL, NULL, NULL, NULL, false);
}

ClusterInstance::ClusterInstance(Forest *f1, Forest *f2, SPRNode *f1_cluster_SPRNode, SPRNode *f2_cluster_SPRNode,bool f2_has_component_zero)
{
    init (f1, f2, f1_cluster_SPRNode, f2_cluster_SPRNode, f2_has_component_zero);
}

ClusterInstance::ClusterInstance(const ClusterInstance &c)
{
    init(c.F1,c.F2,c.F1_cluster_SPRNode,c.F2_cluster_SPRNode,c.F2_has_component_zero);
}

ClusterInstance::~ClusterInstance()
{
}

void ClusterInstance::init(Forest *f1, Forest *f2, SPRNode *f1_cluster_SPRNode,
          SPRNode *f2_cluster_SPRNode, bool f2_has_component_zero)
{
    F1 = f1;
    F2 = f2;
    F1_cluster_SPRNode = f1_cluster_SPRNode;
    F2_cluster_SPRNode = f2_cluster_SPRNode;
    F2_has_component_zero = f2_has_component_zero;
}

bool ClusterInstance::is_original()
{
    if (F1_cluster_SPRNode == NULL && F2_cluster_SPRNode == NULL)
        return true;
    else
        return false;
}

int ClusterInstance::join_cluster(Forest *original_F1, Forest *original_F2)
{
    int adjustment = 0;
    int rho_1 = -1;
    int rho_2 = -1;
    int start = 0;
    Forest *F1_cluster_SPRNode_forest = NULL;
    //Forest *F2_cluster_SPRNode_forest = NULL;
    if (F1_cluster_SPRNode != NULL) {
        F1_cluster_SPRNode->decrease_clustered_children();
        F1_cluster_SPRNode_forest = F1_cluster_SPRNode->get_forest();
    }
    if (F2_cluster_SPRNode != NULL)
        F2_cluster_SPRNode->decrease_clustered_children();
    
    if (F1->contains_rho()) {
        bool contract = true;
        if (F1_cluster_SPRNode->is_leaf()) {
            if (F1_cluster_SPRNode->get_num_clustered_children() >= 1) {
                contract = false;
            }
            else {
                F1_cluster_SPRNode->add_child(F1->get_component(0));
                start = 1;
            }
        }
        else {
            F1_cluster_SPRNode->get_forest()->add_component(F1->get_component(0));
            start = 1;
        }
        if (F1_cluster_SPRNode->parent() == NULL) {
            if (F1_cluster_SPRNode->lchild() != NULL &&
                F1_cluster_SPRNode->lchild()->get_num_clustered_children() > 0) {
                F1_cluster_SPRNode_forest->update_component(
                                                            F1_cluster_SPRNode, F1_cluster_SPRNode->lchild());
            }
            if (F1_cluster_SPRNode->rchild() != NULL &&
                F1_cluster_SPRNode->rchild()->get_num_clustered_children() > 0) {
                F1_cluster_SPRNode_forest->update_component(
                                                            F1_cluster_SPRNode, F1_cluster_SPRNode->rchild());
            }
        }
        if (contract)
            F1_cluster_SPRNode->contract();
    }
    else {
        F1_cluster_SPRNode->add_child(F1->get_component(0));
        if (F1_cluster_SPRNode->get_num_clustered_children() <= 0)
            F1_cluster_SPRNode->contract();
        start = 1;
    }
    // should we add these to a finished_components or something?
    for(int i = start; i < F1->num_components(); i++) {
        if (F1->get_component(i)->str() != "p")
            original_F1->add_component(F1->get_component(i));
        else {
            rho_1 = i;
        }
    }
    
    bool skip = false;
    SPRNode *F2_cluster = F1->get_component(0)->get_twin();
    // TODO: this is still not quite right. We should only add the cluster
    //to the front if that is where it came from
    if ((F1->contains_rho() || F2->contains_rho()) && !F2_has_component_zero && F2_cluster_SPRNode != NULL) {
        bool contract = true;
        if (F2_cluster_SPRNode->is_leaf()) {
            if (F2_cluster_SPRNode->get_num_clustered_children() >= 1) {
                contract = false;
            }
            else {
                //F2_cluster_SPRNode->add_child(F2->get_component(0));
                F2_cluster_SPRNode->add_child(F2_cluster);
                skip = true;
            }
        }
        if (F2_cluster_SPRNode->parent() == NULL) {
            if (F2_cluster_SPRNode->lchild() != NULL &&
                F2_cluster_SPRNode->lchild()->get_num_clustered_children() > 0) {
                F2_cluster_SPRNode->get_forest()->update_component(
                                                                   F2_cluster_SPRNode, F2_cluster_SPRNode->lchild());
            }
            if (F2_cluster_SPRNode->rchild() != NULL &&
                F2_cluster_SPRNode->rchild()->get_num_clustered_children() > 0) {
                F2_cluster_SPRNode->get_forest()->update_component(
                                                                   F2_cluster_SPRNode, F2_cluster_SPRNode->rchild());
            }
        }
        if (contract) {
            F2_cluster_SPRNode->contract();
            if (F1_cluster_SPRNode_forest->get_twin()->get_component(0)->str() == F2_cluster->str() && F1_cluster_SPRNode_forest->get_component(0)->str() == F2_cluster->str() && !F1_cluster_SPRNode_forest->contains_rho()) {
                F1_cluster_SPRNode_forest->add_rho();
                F1_cluster_SPRNode_forest->get_twin()->add_rho();
                adjustment++;
            }
        }
    }
    else { // if (F2_cluster != NULL) {
        if (F2_cluster == NULL)
            F2_cluster = F2->get_component(0);
        skip = true;
        //			cout << "skip=" << skip << endl;
        // TODO: this is not right yet
        if (F2_has_component_zero) {
            //				cout << __LINE__ << endl;
            F1_cluster_SPRNode_forest->get_twin()->add_component(0, F2_cluster);
            if (F1->contains_rho() || F2->contains_rho()) {
                if (!F1_cluster_SPRNode_forest->contains_rho())
                    F1_cluster_SPRNode_forest->add_rho();
                if (!F1_cluster_SPRNode_forest->get_twin()->contains_rho())
                    F1_cluster_SPRNode_forest->get_twin()->add_rho();
            }
        }
        else if (F2_cluster_SPRNode == NULL) {
            skip = false;
        }
        else {
            // TODO: check this
            F2_cluster_SPRNode->add_child(F2->get_component(0));
            F2_cluster = F2->get_component(0);
            if ((F2_cluster_SPRNode->lchild() == NULL || F2_cluster_SPRNode->rchild() == NULL) && F2_cluster_SPRNode->get_num_clustered_children() <= 0) {
                F2_cluster_SPRNode->contract();
            }
        }
    }
    // should we add these to a finished_components or something?
    for(int i = 0; i < F2->num_components(); i++) {
        if (F2->get_component(i) == F2_cluster) {
            if (!skip) {
                if (F2->get_component(i)->str() != "p")
                    F1_cluster_SPRNode_forest->get_twin()->add_component(F2->get_component(i));
                else
                    rho_2 = i;
            }
        }
        else {
            if (F2->get_component(i)->str() != "p")
                original_F2->add_component(F2->get_component(i));
            else
                rho_2 = i;
        }
    }
    
    if (rho_1 >= 0) {
        F1->get_component(rho_1)->delete_tree();
        F1->set_component(rho_1,NULL);
    }
    if (rho_2 >= 0) {
        F2->get_component(rho_2)->delete_tree();
        F2->set_component(rho_2,NULL);
    }
    F1->erase_components();
    F2->erase_components();
    return adjustment;
}

void cluster_reduction_find_components(SPRNode *n,
                                       std::vector<bool> *F2_cluster_copy_components,
                                       std::vector<bool> *old_F2_keep_components,
                                       int cluster_component_number)
{
    SPRNode *lc = n->lchild();
    SPRNode *rc = n->rchild();
    if (lc != NULL)
        cluster_reduction_find_components(lc, F2_cluster_copy_components,
                                          old_F2_keep_components, cluster_component_number);
    if (rc != NULL)
        cluster_reduction_find_components(rc, F2_cluster_copy_components,
                                          old_F2_keep_components, cluster_component_number);
    if (lc == NULL && rc == NULL) {
        int cnumber = n->get_twin()->get_component_number();
        if (cnumber != cluster_component_number) {
            (*F2_cluster_copy_components)[cnumber] = true;
            (*old_F2_keep_components)[cnumber] = false;
        }
    }
}

std::list<ClusterInstance> cluster_reduction(Forest *old_F1, Forest *old_F2,
                                             std::list<SPRNode *> *cluster_points)
{
    std::list<ClusterInstance> clusters = std::list<ClusterInstance>();
    std::vector<bool> old_F2_keep_components =
    std::vector<bool>(old_F2->num_components(), true);
    for(std::list<SPRNode *>::iterator i = cluster_points->begin();
        i != cluster_points->end(); i++) {
        // Cluster F1
        SPRNode *F1_root_SPRNode = *i;
        SPRNode *F1_cluster_SPRNode = F1_root_SPRNode->parent();
        
        F1_root_SPRNode->cut_parent();
        std::vector<SPRNode *> cluster_F1_components = std::vector<SPRNode *>();
        cluster_F1_components.push_back(F1_root_SPRNode);
        Forest *F1 = new Forest(cluster_F1_components);
        
        // Cluster F2
        SPRNode *F2_root_SPRNode = F1_root_SPRNode->get_twin();
        SPRNode *F2_cluster_SPRNode = F2_root_SPRNode->parent();
        std::vector<bool> F2_cluster_copy_components =
        std::vector<bool>(old_F2->num_components(), false);
        int cnumber = F2_root_SPRNode->get_component_number();
        bool F2_has_component_zero = false;
        bool skip_F2_cluster = false;
        if (F2_root_SPRNode->parent() != NULL) {
            F2_root_SPRNode->cut_parent();
        }
        else {
            if (old_F2_keep_components[cnumber] == false) {
                skip_F2_cluster = true;
            }
            else if (old_F2->get_component(cnumber) == F2_root_SPRNode){
                old_F2_keep_components[cnumber] = false;
                if (cnumber == 0) {
                    F2_has_component_zero = true;
                }
            }
            else {
                //				if (F2_root_SPRNode->get_forest()->get_cluster()->F2_has_component_zero == false) {
                F2_root_SPRNode->get_forest()->get_cluster()
                ->F2_has_component_zero = true;
                //					if (F2_root_SPRNode->get_forest()->contains_rho() == false)
                //						F2_root_SPRNode->get_forest()->add_rho();
                F2_cluster_SPRNode = F2_root_SPRNode->get_forest()->get_cluster()->F2_cluster_SPRNode;
                //				}
                skip_F2_cluster = true;
            }
        }
        
        cluster_reduction_find_components(F1_root_SPRNode,
                                          &F2_cluster_copy_components, &old_F2_keep_components, cnumber);
        std::vector<SPRNode *> cluster_F2_components = std::vector<SPRNode *>();
        if (!skip_F2_cluster)
            cluster_F2_components.push_back(F2_root_SPRNode);
        for(int i = 0; i < F2_cluster_copy_components.size(); i++) {
            if (F2_cluster_copy_components[i] == true)
                cluster_F2_components.push_back(old_F2->get_component(i));
        }
        Forest *F2 = new Forest(cluster_F2_components);
        clusters.push_back(ClusterInstance(F1, F2, F1_cluster_SPRNode,
                                           F2_cluster_SPRNode, F2_has_component_zero));
        F1->set_cluster(clusters.back());
        F2->set_cluster(clusters.back());
        if (F1_cluster_SPRNode != NULL)
            F1_cluster_SPRNode->increase_clustered_children();
        if (F2_cluster_SPRNode != NULL)
            F2_cluster_SPRNode->increase_clustered_children();
        F1->label_SPRNodes_with_forest();
        F2->label_SPRNodes_with_forest();
        F1->set_twin(F2);
        F2->set_twin(F1);
        
        for(int i = 0; i < F2_cluster_copy_components.size(); i++) {
            if (F2_cluster_copy_components[i] == true) {
                cluster_F2_components.push_back(
                                                (old_F2->get_component(i)));
            }
        }
    }
    // remove any clustered components from old_F2
    std::vector<SPRNode *> old_F2_remaining_components = std::vector<SPRNode *>();
    for(int i = 0; i < old_F2_keep_components.size(); i++) {
        if (old_F2_keep_components[i] == true) {
            old_F2_remaining_components.push_back(
                                                  (old_F2->get_component(i)));
            }
    }
    Forest *replace_old_F2 = new Forest(old_F2_remaining_components);
    replace_old_F2->swap(old_F2);
    replace_old_F2->erase_components();
    delete replace_old_F2;
    
    
    clusters.push_back(ClusterInstance(old_F1, old_F2, NULL, NULL, true));
    old_F1->set_cluster(clusters.back());
    old_F2->set_cluster(clusters.back());
    old_F1->label_SPRNodes_with_forest();
    old_F2->label_SPRNodes_with_forest();
    old_F1->set_twin(old_F2);
    old_F2->set_twin(old_F1);
    
    return clusters;
    
}


