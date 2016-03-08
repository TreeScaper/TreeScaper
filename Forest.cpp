#include "Forest.hpp"

Forest::Forest()
{
    init(std::vector<SPRNode *>());
}

Forest::Forest(std::vector<SPRNode *> components)
{
    init(components);
}

Forest::Forest(SPRNode *head)
{
    components = std::vector<SPRNode *>();
    components.push_back(new SPRNode(*head));
    deleted_nodes = std::vector<SPRNode *>();
    rho = false;
    twin = NULL;
    cluster = NULL;
}

Forest::Forest(const Forest &f)
{
    components = std::vector<SPRNode *>(f.components.size());
    for(int i = 0; i < f.components.size(); i++) {
        //if (f.components[i] != NULL)
        components[i] = new SPRNode(*f.components[i]);
    }
    deleted_nodes = std::vector<SPRNode *>();
    rho = f.rho;
    twin = NULL;
    cluster = f.cluster;
    //label_SPRNodes_with_forest();
}

Forest::Forest(Forest *f)
{
    components = std::vector<SPRNode *>(f->components.size());
    for(int i = 0; i < f->components.size(); i++) {
        //if (f->components[i] != NULL)
        components[i] = new SPRNode(*f->components[i]);
    }
    deleted_nodes = std::vector<SPRNode *>();
    rho = f->rho;
    twin = NULL;
    cluster = f->cluster;
    //label_SPRNodes_with_forest();
}

Forest::Forest(Forest *f, bool b)
{
    components = std::vector<SPRNode *>(f->components.size());
    for(int i = 0; i < f->components.size(); i++) {
        //if (f->components[i] != NULL)
        components[i] = new SPRNode(*f->components[i]);
    }
    deleted_nodes = std::vector<SPRNode *>();
    rho = f->rho;
    twin = NULL;
    cluster = f->cluster;
    //label_SPRNodes_with_forest();
}


void Forest::init(std::vector<SPRNode *> components)
{
    this->components = std::vector<SPRNode *>(components);
    deleted_nodes = std::vector<SPRNode *>();
    rho = false;
    for(int i = 0; i < components.size(); i++) {
        if (components[i]->str() == "p")
            rho = true;
    }
    twin = NULL;
    cluster = NULL;
}

Forest::~Forest()
{
    for(int i = 0; i < components.size(); i++) {
        //if (components[i] != NULL) {
        components[i]->delete_tree();
        components[i] = NULL;
        //}
    }
    for(int i = 0; i < deleted_nodes.size(); i++) {
        //if (deleted_SPRNodes[i] != NULL) {
        deleted_nodes[i]->delete_tree();
        deleted_nodes[i] = NULL;
        //}
    }
}

// swap the contents of two forests
void Forest::swap(Forest *f)
{
    std::vector<SPRNode *> components_temp = this->components;
    this->components = f->components;
    f->components = components_temp;

    std::vector<SPRNode *> deleted_nodes_temp = this->deleted_nodes;
    this->deleted_nodes = f->deleted_nodes;
    f->deleted_nodes = deleted_nodes_temp;
    
    bool rho_temp = this->rho;
    this->rho = f->rho;
    f->rho = rho_temp;
    
    ClusterInstance *c_temp = this->cluster;
    this->cluster = f->cluster;
    f->cluster = c_temp;
}

// print the forest
void Forest::print_components()
{
    std::vector<SPRNode *>::iterator it = components.begin();
    for(it = components.begin(); it != components.end(); it++) {
        SPRNode *root = *it;
        if (root == NULL)
            std::cout << "!";
        else if (root->is_leaf() && root->str() == "")
            std::cout << "*";
        else
            root->print_subtree_hlpr();
        std::cout << " ";
    }
    std::cout << std::endl;
}

// print the components seperated by s
void Forest::print_components(std::string s)
{
    std::vector<SPRNode *>::iterator it = components.begin();
    for(it = components.begin(); it != components.end(); it++) {
        SPRNode *root = *it;
        if (root == NULL)
            std::cout << "!";
        else
            root->print_subtree_hlpr();
        std::cout << s;
    }
    std::cout << std::endl;
}

// print the forest showing twins
void Forest::print_components_with_twins()
{
    std::vector<SPRNode *>::iterator it = components.begin();
    for(it = components.begin(); it != components.end(); it++) {
        SPRNode *root = *it;
        if (root == NULL)
            std::cout << "!";
        else
            root->print_subtree_twin_hlpr();
        std::cout << " ";
    }
    std::cout << std::endl;
}

// print the forest showing edge intervals
void Forest::print_components_with_edge_pre_interval()
{
    std::vector<SPRNode *>::iterator it = components.begin();
    for(it = components.begin(); it != components.end(); it++) {
        SPRNode *root = *it;
        if (root == NULL)
            std::cout << "!";
        else
            std::cout << root->str_edge_pre_interval_subtree();
        std::cout << " ";
    }
    std::cout << std::endl;
}

// return the string for this forest
std::string Forest::str()
{
    std::string s = "";
    std::vector<SPRNode *>::iterator it;
    for(it = components.begin(); it != components.end(); it++) {
        SPRNode *root = *it;
        if (root == NULL)
            s += "!";
        else
            s += root->str_subtree();
        s += " ";
    }
    return s;
}

void Forest::add_component(SPRNode *head)
{
    components.push_back(head);
}

void Forest::add_component(int pos, SPRNode *head)
{
    components.insert(components.begin() + pos, head);
}

void Forest::add_deleted_node(SPRNode *n)
{
    deleted_nodes.push_back(n);
}

SPRNode *Forest::get_component(int i)
{
    return components[i];
}

void Forest::set_component(int i, SPRNode *head)
{
    components[i] = head;
}

inline size_t Forest::num_components()
{
    return components.size();
}

void Forest::set_twin(Forest *f)
{
    twin = f;
}

Forest *Forest::get_twin()
{
    return twin;
}

ClusterInstance *Forest::get_cluster()
{
    return cluster;
}

void Forest::set_cluster(ClusterInstance *c)
{
    cluster = c;
}

void Forest::set_cluster(ClusterInstance &c)
{
    cluster = &c;
}

void Forest::update_component(SPRNode *old_c, SPRNode *new_c)
{
    std::vector<SPRNode *>::iterator i;
    for(i = components.begin(); i != components.end(); i++) {
        SPRNode *component = *i;
        if (&(*component) == &(*old_c))
            *i = new_c;
    }
}

// return a list of the sibling pairs
std::list<SPRNode *> *Forest::find_sibling_pairs()
{
    std::list<SPRNode *> *sibling_pairs = new std::list<SPRNode *>();
    std::vector<SPRNode *>::iterator i;
    for(i = components.begin(); i != components.end(); i++) {
        SPRNode *component = *i;
        component->append_sibling_pairs(sibling_pairs);
    }
    return sibling_pairs;
}

// return a deque of the singleton leaves
std::list<SPRNode *> Forest::find_singletons()
{
    std::list<SPRNode *> singletons = std::list<SPRNode *>();
    std::vector<SPRNode *>::iterator i = components.begin();
    // TODO: is this a hack or correct? We don't want the first
    // component to be a singleton because of rho!
    i++;
    for(; i != components.end(); i++) {
        SPRNode *component = *i;
        if (component->is_leaf()) {
            singletons.push_back(component);
        }
    }
    return singletons;
}

// make SPRNodes pointed to in the forest point back
void Forest::resync()
{
    std::vector<SPRNode *>::iterator i;
    for(i = components.begin(); i != components.end(); i++) {
        (*i)->resync();
    }
}

// clear twin pointers
void Forest::unsync()
{
    std::vector<SPRNode *>::iterator i;
    for(i = components.begin(); i != components.end(); i++) {
        (*i)->unsync();
    }
}

// clear interior twin pointers
void Forest::unsync_interior()
{
    std::vector<SPRNode *>::iterator i;
    for(i = components.begin(); i != components.end(); i++) {
        (*i)->unsync_interior();
    }
}

// clear twin pointers
SPRNode *Forest::find_by_prenum(int prenum)
{
    std::vector<SPRNode *>::iterator i;
    for(i = components.begin(); i != components.end(); i++) {
        SPRNode *ans = (*i)->find_by_prenum(prenum);
        if (ans != NULL)
            return ans;
    }
    return NULL;
}

void Forest::labels_to_numbers(std::map<std::string, int> *label_map, std::map<int, std::string> *reverse_label_map)
{
    std::vector<SPRNode *>::iterator i;
    for(i = components.begin(); i != components.end(); i++)
    {
        (*i)->labels_to_numbers(label_map, reverse_label_map);
    }
}

void Forest::numbers_to_labels(std::map<int, std::string> *reverse_label_map)
{
    std::vector<SPRNode *>::iterator i;
    for(i = components.begin(); i != components.end(); i++)
    {
        (*i)->numbers_to_labels(reverse_label_map);
    }
}

size_t Forest::size()
{
    return components.size();
}

bool Forest::add_rho()
{
    if (rho)
        return false;
    SPRNode *T_p = new SPRNode("p");
    add_component(T_p);
    rho = true;
    return true;
}

void Forest::set_rho(bool b)
{
    rho = b;
}

bool Forest::contains_rho()
{
    return rho;
}

void Forest::erase_components(int start, int end)
{
    components.erase(components.begin()+start, components.begin()+end);
}

void Forest::erase_components()
{
    components.clear();
}

// tell the SPRNodes what forest they are in
void Forest::label_SPRNodes_with_forest()
{
    for(int i = 0; i < size(); i++) {
        get_component(i)->set_forest_rec(this);
    }
}

void Forest::move_first_component_to_end()
{
    SPRNode *temp = components[0];
    components[0] = components[components.size()-1];
    components[components.size()-1] = temp;
}

void Forest::unprotect_edges()
{
    for(int i = 0; i < size(); i++)
    {
        get_component(i)->unprotect_subtree();
    }
}

// Make the leaves of two forests point to their twin in the other tree
// Note: removes unique leaves
bool sync_twins(Forest *T1, Forest *T2)
{
    std::vector<SPRNode *> T1_labels = std::vector<SPRNode *>();
    std::vector<SPRNode *> T2_labels = std::vector<SPRNode *>();
    std::vector<SPRNode *> T1_components = T1->components;
    std::vector<SPRNode *> T2_components = T2->components;
    std::vector<SPRNode *>::iterator i;
    SPRNode *T1_rho = NULL;
    SPRNode *T2_rho = NULL;
    for(i = T1_components.begin(); i != T1_components.end(); i++) {
        SPRNode *component = *i;
        std::vector<SPRNode *> unsorted_labels = component->find_leaves();
        std::vector<SPRNode *>::iterator j;
        for(j = unsorted_labels.begin(); j != unsorted_labels.end(); j++) {
            SPRNode *leaf = *j;
            //			cout << "T1: " << leaf->str() << endl;
            if (leaf->str() == "p") {
                T1_rho = leaf;
            }
            else {
                // find smallest number contained in the label
                int number = stomini(leaf->str());
                //				cout << "\t" << number << endl;
                if (number < INT_MAX) {
                    if (number >= T1_labels.size())
                        T1_labels.resize(number+1, 0);
                    T1_labels[number] = leaf;
                }
            }
        }
    }
    for(i = T2_components.begin(); i != T2_components.end(); i++) {
        SPRNode *component = *i;
        std::vector<SPRNode *> unsorted_labels = component->find_leaves();
        std::vector<SPRNode *>::iterator j;
        for(j = unsorted_labels.begin(); j != unsorted_labels.end(); j++) {
            SPRNode *leaf = *j;
            if (leaf->str() == "p") {
                T2_rho = leaf;
            }
            else {
                // find smallest number contained in the label
                int number = stomini(leaf->str());
                if (number < INT_MAX) {
                    if (number >= T2_labels.size())
                        T2_labels.resize(number+1, 0);
                    T2_labels[number] = leaf;
                }
            }
        }
    }
    T1_labels.resize(T1_labels.size()+1);
    T1_labels[T1_labels.size()-1]=T1_rho;
    T2_labels.resize(T2_labels.size()+1);
    T2_labels[T2_labels.size()-1]=T2_rho;
    
    size_t size = T1_labels.size();
    if (size > T2_labels.size())
        size = T2_labels.size();
    //	cout << "Syncing Twins" << endl;
    for(int i = 0; i < size; i++) {
        SPRNode *T1_a = T1_labels[i];
        SPRNode *T2_a = T2_labels[i];
        if (T1_a == NULL && T2_a != NULL) {
            SPRNode *node = T2_a->parent();
            if (node == NULL)
                return false;
            size_t numc = node->get_children().size();
            if (node->parent() == NULL && node->lchild()->is_leaf() &&
                (numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
                return false;
                SPRNode *sibling = node->lchild();
                if (sibling == T2_a)
                    sibling = node->rchild();
                T2_labels[stomini(sibling->str())] = sibling;
            }
            delete T2_a;
            if (node->get_children().size() < 2) {
                if (node->get_children().size() == 1)
                    node->lchild()->lost_child();
                node = node->contract(true);
            }
        }
        else if (T2_a == NULL && T1_a != NULL) {
            SPRNode *node = T1_a->parent();
            if (node == NULL)
                return false;
            size_t numc = node->get_children().size();
            if (node->parent() == NULL && node->lchild()->is_leaf() &&
                (numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
                return false;
                SPRNode *sibling = node->lchild();
                if (sibling == T1_a)
                    sibling = node->rchild();
                T1_labels[stomini(sibling->str())] = sibling;
            }
            delete T1_a;
            if (node->get_children().size() < 2) {
                if (node->get_children().size() == 1)
                    node->lchild()->lost_child();
                node = node->contract(true);
            }
            
        }
        if (T1_a != NULL && T2_a != NULL) {
            T1_a->set_twin(T2_a);
            T2_a->set_twin(T1_a);
            //			cout << T1_a->str() << endl;
        }
    }
    for(size_t i = size; i < T1_labels.size(); i++) {
        SPRNode *T1_a = T1_labels[i];
        if (T1_a != NULL) {
            SPRNode *node = T1_a->parent();
            if (node == NULL)
                return false;
            size_t numc = node->get_children().size();
            if (node->parent() == NULL && node->lchild()->is_leaf() &&
                (numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
                return false;
                SPRNode *sibling = node->lchild();
                if (sibling == T1_a)
                    sibling = node->rchild();
                T1_labels[stomini(sibling->str())] = sibling;
            }
            delete T1_a;
            if (node->get_children().size() < 2) {
                if (node->get_children().size() == 1)
                    node->lchild()->lost_child();
                node = node->contract(true);
            }
            
        }
    }
    for(size_t i = size; i < T2_labels.size(); i++) {
        SPRNode *T2_a = T2_labels[i];
        if (T2_a != NULL) {
            SPRNode *node = T2_a->parent();
            if (node == NULL)
                return false;
            size_t numc = node->get_children().size();
            if (node->parent() == NULL && node->lchild()->is_leaf() &&
                (numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
                return false;
                SPRNode *sibling = node->lchild();
                if (sibling == T2_a)
                    sibling = node->rchild();
                T2_labels[stomini(sibling->str())] = sibling;
            }
            delete T2_a;
            if (node->get_children().size() < 2) {
                if (node->get_children().size() == 1)
                    node->lchild()->lost_child();
                node = node->contract(true);
            }
        }
    }
    return true;
}

/* make interior SPRNodes point to the lca of their descendants in the other tree
 assumes that sync_twins has already been called
 assumes that component 1 of T1 matches with 1 of T2
 NOTE: this isn't true during the algorithm so this will need to be changed
 if we want to interleave clustering. It should be just component 1 of T1
 matching multiple components of T2 (The first several components?)
 */
void sync_interior_twins(Forest *T1, Forest *T2)
{
    SPRNode  *root1 = T1->get_component(0);
    SPRNode  *root2 = T2->get_component(0);
    LCA T1_LCA = LCA(root1);
    LCA T2_LCA = LCA(root2);
    sync_interior_twins(root1, &T2_LCA);
    sync_interior_twins(root2, &T1_LCA);
}

/* make interior SPRNodes point to the LCA of their descendants in the other
 forest if there is one unambiguous LCA
 * This is true for a SPRNode n of T1 if all leaves that are a descendant
 of T1 either map to a single component of F2 or are from another
 component of F2 such that the root of that component maps to
 a descendant of n (i.e. a finished component)
 * assumes that sync_twins has already been called
 */
// TODO: initializing parameters seems to be slow
void sync_interior_twins_real(Forest *T1, Forest *F2)
{
    SPRNode  *T1_root = T1->get_component(0);
    LCA T1_LCA = LCA(T1_root);
    //int T1_size = T1_root->size_using_prenum();
    // roots of F2
    std::vector<SPRNode *> F2_roots = std::vector<SPRNode *>();
    // LCA queries for F2
    std::vector<LCA> F2_LCAs = std::vector<LCA>();
    // lists of root SPRNodes that map to a given T1 SPRNode
    T1_root->initialize_root_lcas(std::list<SPRNode *>());
    // list of active descendants
    T1_root->initialize_active_descendants(std::list<SPRNode *>());
    
    // should be fine.
    for(int i = 0; i < F2->num_components(); i++) {
        //		cout << "starting i" << endl;
        F2_roots.push_back(F2->get_component(i));
        if (F2_roots[i]->str() == "p") {
            //	F2_LCAs.push_back(LCA());
            F2_LCAs.push_back(F2_roots[i]);
            //F2_LCAs.push_back(NULL);//LCA(F2_roots[i]));
            continue;
        }
        F2_LCAs.push_back(LCA(F2_roots[i]));
        // number the component
        F2_roots[i]->initialize_component_number(i);
        // list of SPRNodes that get deleted when a component is finished
        F2_roots[i]->initialize_removable_descendants(std::list<std::list<SPRNode *>::iterator>());
        // sync the component with T1
        if (F2_roots[i]->str() != "p" &&
            !(F2_roots[i]->get_twin() != NULL && F2_roots[i]->get_twin()->parent() == NULL)) {
            sync_interior_twins(F2_roots[i], &T1_LCA);
        }
        if (i > 0 || T1->contains_rho())
            F2_roots[i]->get_twin()->get_root_lcas()->push_back(F2_roots[i]);
        else
            T1_root->get_root_lcas()->push_back(F2_roots[i]);
        //		cout << "b" << endl;
    }
    sync_interior_twins(T1_root, &F2_LCAs);
}

/* make interior SPRNodes point to the lca of their descendants in the other
 * tree
 * assumes that sync_twins has already been called
 */
void sync_interior_twins(SPRNode *n, LCA *twin_LCA)
{
    std::list<SPRNode *>::iterator c = n->get_children().begin();
    if (c == n->get_children().end())
        return;
    sync_interior_twins(*c, twin_LCA);
    n->set_twin((*c)->get_twin());
    c++;
    while(c != n->get_children().end()) {
        sync_interior_twins(*c, twin_LCA);
        SPRNode *twin = twin_LCA->get_lca(n->get_twin(), (*c)->get_twin());
        n->set_twin(twin);
        c++;
    }
}

void sync_interior_twins(SPRNode *n, std::vector<LCA> *F2_LCAs)
{
    SPRNode *lc = n->lchild();
    SPRNode *rc = n->rchild();
    std::list<SPRNode *> *active_descendants = n->get_active_descendants();
    // visit children first
    std::list<SPRNode *>::iterator c;
    for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
        sync_interior_twins(*c, F2_LCAs);
    }
    if (n->get_children().size() == 0) {
        active_descendants->push_back(n->get_twin());
        std::list<SPRNode *>::iterator SPRNode_location = active_descendants->end();
        SPRNode_location--;
        n->get_twin()->get_removable_descendants()->push_back(SPRNode_location);
    }
    // no rc so propogate up
    if (n->get_children().size() == 1) {
        SPRNode *lc = n->get_children().front();
        n->set_twin(lc->get_twin());
        std::list<SPRNode *> *lc_active_descendants = lc->get_active_descendants();
        active_descendants->splice(active_descendants->end(),*lc_active_descendants);
    }
    // TODO: generalize from here for 2 or more children
    // two children so put their info together
    else if (lc != NULL && rc != NULL) {
        std::vector<std::list<SPRNode *>::iterator> SPRNode_location =
        std::vector<std::list<SPRNode *>::iterator>();
        std::list<SPRNode *>::iterator SPRNode1_location;
        int nonempty_active_descendants_count = 0;
        for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
            if (!(*c)->get_active_descendants()->empty()) {
                nonempty_active_descendants_count++;
                if (nonempty_active_descendants_count > 1) {
                    SPRNode1_location = active_descendants->end();
                    SPRNode1_location--;
                    SPRNode_location.push_back(SPRNode1_location);
                }
                //		cout << active_descendants->size() << endl;
                active_descendants->splice(active_descendants->end(),
                                           *((*c)->get_active_descendants()));
                //		cout << active_descendants->size() << endl;
            }
        }
        
        /* check the intersection points to see if we have two
         leaves from the same component
         */
        for(int i = 0; i < SPRNode_location.size(); i++) {
            std::list<SPRNode *>::iterator SPRNode1_location = SPRNode_location[i];
            std::list<SPRNode *>::iterator SPRNode2_location = SPRNode1_location;
            SPRNode2_location++;
            delete_and_merge_LCAs(active_descendants, F2_LCAs, SPRNode1_location,
                                  SPRNode2_location);
        }
        
        /* check to see if n is twinned by a root of F2
         if so, then remove each leaf twinned by that component
         and check each of the new intersection points
         */
        std::list<SPRNode *> *root_lcas = n->get_root_lcas();
        while(!root_lcas->empty()) {
            
            SPRNode *root_lca = root_lcas->front();
            root_lcas->pop_front();
            /* TODO: problem when n is a root
             We don't care about this but it might mean there is a different
             problem
             */
            if (n->parent() != NULL) {
                delete_and_merge_LCAs(root_lca, active_descendants, F2_LCAs);
            }
        }
        /* If we have a single element in n's active descendants
         list then set twin pointers appropriately
         */
    
        if (active_descendants->size() == 1) {
            SPRNode *twin = active_descendants->front();
            n->set_twin(twin);
        }
    }
    if (n->parent() == NULL)
        active_descendants->clear();
}

void sync_af_twins(Forest *F1, Forest *F2)
{
    F1->unsync();
    F2->unsync();
    sync_twins(F1, F2);
    for(int i = 0; i < F1->num_components(); i++)
    {
        F1->get_component(i)->sync_af_twins();
    }
}

/* merge two SPRNodes from a list into their LCA if they are from
 the same component
 */
void delete_and_merge_LCAs(std::list<SPRNode *> *active_descendants,
                           std::vector<LCA> *F2_LCAs, std::list<SPRNode *>:: iterator SPRNode1_location,
                           std::list<SPRNode *>:: iterator SPRNode2_location)
{
    //	while(active_descendants->size() > 1) {
    SPRNode *n1 = *SPRNode1_location;
    SPRNode *n2 = *SPRNode2_location;
    int component1 = n1->get_component_number();
    int component2 = n2->get_component_number();
    if (component1 == component2) {
        SPRNode *lca = (*F2_LCAs)[component1].get_lca(n1,n2);
        
        //		cout << "xa" << endl;
        std::list<SPRNode *>::iterator lca_location =
        active_descendants->insert(SPRNode1_location,lca);
        //		cout << "xb" << endl;
        //		active_descendants->erase(SPRNode1_location);
        //		cout << "xc" << endl;
        
        // TODO: could this be faster?
        bool remove = false;
        std::list<std::list<SPRNode *>::iterator>::iterator i;
        for(i = n1->get_removable_descendants()->begin(); i != n1->get_removable_descendants()->end(); i++) {
            if (*i == SPRNode1_location) {
                //active_descendants->erase(*i);
                remove = true;
                break;
            }
        }
        if (remove) {
            active_descendants->erase(*i);
            n1->get_removable_descendants()->erase(i);
        }
        //		n1->get_removable_descendants()->clear();
        remove = false;
        for(i = n2->get_removable_descendants()->begin(); i != n2->get_removable_descendants()->end(); i++) {
            if (*i == SPRNode2_location) {
                //active_descendants->erase(*i);
                remove = true;
                break;
            }
        }
        if (remove) {
            active_descendants->erase(*i);
            n2->get_removable_descendants()->erase(i);
        }
        lca->get_removable_descendants()->push_back(lca_location);
    }
}

/* delete each leaf from the list that is twinned with the component
 of n. For each such deleted SPRNode, merge its predecessor
 and successor in the list into their LCA if they are from
 the same component (other than n's component)
 */
void delete_and_merge_LCAs(SPRNode *n, std::list<SPRNode *>
                           *active_descendants, std::vector<LCA> *F2_LCAs)
{
    int component = n->get_component_number();
    std::list<std::list<SPRNode *>::iterator> *removable_descendants	=
    n->get_removable_descendants();
    
    if (n->lchild() != NULL)
        delete_and_merge_LCAs(n->lchild(), active_descendants, F2_LCAs);
    if (n->rchild() != NULL)
        delete_and_merge_LCAs(n->rchild(), active_descendants, F2_LCAs);

    while (!removable_descendants->empty())
    {
        std::list<SPRNode *>::iterator leaf_location = removable_descendants->front();
        removable_descendants->pop_front();
        
        // TODO: problem here
        if (leaf_location != active_descendants->begin() &&
            leaf_location != active_descendants->end() &&
            leaf_location != -- active_descendants->end()) {
            //		if (active_descendants->front() != *leaf_location
            //				&& active_descendants->back() != *leaf_location) {
            std::list<SPRNode *>::iterator SPRNode1_location = leaf_location;
            //		cout << "fooe" << endl;
            std::list<SPRNode *>::iterator SPRNode2_location = leaf_location;
            //		cout << "foof" << endl;
            SPRNode1_location--;
            //		cout << "foog" << endl;
            SPRNode2_location++;
            //		cout << "fooh" << endl;
            active_descendants->erase(leaf_location);
            //		cout << "fooi" << endl;
            int SPRNode1_component = (*SPRNode1_location)->get_component_number();
            //		cout << "fooj" << endl;
            if (component != SPRNode1_component)
                delete_and_merge_LCAs(active_descendants, F2_LCAs, SPRNode1_location,
                                      SPRNode2_location);
            //		cout << "fook" << endl;
        }
        else {//if (active_descendants->size() > 1){
            active_descendants->erase(leaf_location);
        }
        
    }
}

std::list<SPRNode *> *find_cluster_points(Forest *F1, Forest *F2)
{
    std::list<SPRNode *> *cluster_points = new std::list<SPRNode *>();
    std::vector<int> *leaf_counts_F1 = NULL;
    std::vector<int> *leaf_counts_F2 = NULL;
    //if (MULTI_CLUSTER) {
    if (false) {
        leaf_counts_F1 = F1->get_component(0)->find_leaf_counts();
        leaf_counts_F2 = F2->get_component(0)->find_leaf_counts();
    }
    find_cluster_points(F1->get_component(0), cluster_points, leaf_counts_F1,
                        leaf_counts_F2);
    //if (MULTI_CLUSTER) {
    if (false) {
        delete leaf_counts_F1;
        //		delete leaf_counts_F2;
    }
    //cout << "foo" << endl;
    return cluster_points;
}

// find the cluster points
void find_cluster_points(SPRNode *n, std::list<SPRNode *> *cluster_points,
                         std::vector<int> *leaf_counts_F1, std::vector<int> *leaf_counts_F2)
{
    std::list<SPRNode *>::iterator c;
    for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
        find_cluster_points(*c, cluster_points, leaf_counts_F1,
                            leaf_counts_F2);
    }
    bool is_cluster = true;
    int num_clustered_children = 0;
    if (n->get_twin() == NULL ||
        n->parent() == NULL ||
        n->get_children().size() < 2 ||
        n->get_SPR_depth() > n->get_twin()->get_twin()->get_SPR_depth())
        is_cluster = false;
    else {
            for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
                if ((*c)->get_twin() != NULL &&
                    (*c)->get_SPR_depth() <= (*c)->get_twin()->get_twin()->get_SPR_depth())
                    num_clustered_children++;
            }
            if (num_clustered_children == n->get_children().size())
                is_cluster = false;
    }
    if (is_cluster) {
        //		cout << "added cluster_point" << endl;
        cluster_points->push_back(n);
    }
    // buggy, needs testing, doesn't seem worth it
    //else if (MULTI_CLUSTER && n->get_twin() != NULL && n->parent() != NULL &&
    else if (false && n->get_twin() != NULL && n->parent() != NULL &&
             n->get_children().size() >= 2) {
        // TODO: use find_leaf_counts if this works
        SPRNode *n_twin = n->get_twin();
        int num_leaves = (*leaf_counts_F1)[n->get_preorder_number()];
        std::vector<SPRNode *> chosen = std::vector<SPRNode *>();
        int chosen_leaves = 0;
        if (n_twin != NULL && n->get_edge_pre_start() > -1 && n->get_edge_pre_end() > -1 && n_twin->get_children().size() > 2) {
            for(c = n_twin->get_children().begin(); c != n_twin->get_children().end(); c++) {
                int c_num_leaves = (*leaf_counts_F2)[(*c)->get_preorder_number()];
                int c_twin_pre = (*c)->get_twin()->get_preorder_number();
                if (c_twin_pre >= n->get_edge_pre_start() &&
                    c_twin_pre <= n->get_edge_pre_end()) {
                    chosen.push_back(*c);
                    chosen_leaves += c_num_leaves;
                }
            }
            // PROBLEM: the new SPRNode should have its own preorder number
            // and its own size
            if (num_leaves == chosen_leaves) {
                SPRNode *new_child = new SPRNode();
                n_twin->add_child(new_child);
                new_child->set_preorder_number(n_twin->get_preorder_number());
                new_child->set_edge_pre_start(n_twin->get_edge_pre_start());
                new_child->set_edge_pre_end(n_twin->get_edge_pre_end());
                n->set_twin(new_child);
                new_child->set_twin(n);
                cluster_points->push_back(n);
                for(int i = 0; i < chosen.size(); i++) {
                    new_child->add_child(chosen[i]);
                }
            }
            
        }
    }
}

// swap two forests
void swap(Forest **a, Forest **b)
{
    (*a)->swap(*b);
}

// expand all contracted SPRNodes
void expand_contracted_SPRNodes(Forest *F)
{
    for(int i = 0; i < F->num_components(); i++) {
        F->get_component(i)->expand_contracted_SPRNodes();
    }
}

Forest *build_finished_forest(std::string &name)
{
    Forest *new_forest = new Forest();
    size_t old_loc = 0;
    size_t loc = 0;
    while ((loc = name.find(" ", old_loc)) != std::string::npos)
    {
        new_forest->add_component(new SPRNode(name.substr(old_loc,loc-old_loc)));
        if (name.substr(old_loc,loc-old_loc) == "p")
            new_forest->set_rho(true);
        //new_forest->print_components();
        old_loc = loc+1;
    }
    new_forest->add_component(new SPRNode(name.substr(old_loc,loc-old_loc)));
    return new_forest;
    //new_forest->add_component(spr_building_tree(name.substr(old_loc,loc-old_loc)));
}

Forest *build_forest(std::string &name)
{
    Forest *new_forest = new Forest();
    size_t old_loc = 0;
    size_t loc = 0;
    while ((loc = name.find(" ", old_loc)) != std::string::npos)
    {
        new_forest->add_component(spr_building_tree(name.substr(old_loc,loc-old_loc)));
        if (name.substr(old_loc,loc-old_loc) == "p")
            new_forest->set_rho(true);
        //new_forest->print_components();
        old_loc = loc+1;
    }
    new_forest->add_component(spr_building_tree(name.substr(old_loc,loc-old_loc)));
    return new_forest;
}
