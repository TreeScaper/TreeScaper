#include "SPRNode.hpp"

//#define COPY_CONTRACTED


SPRNode::SPRNode()
{
    init(NULL, NULL, NULL, "", 0);
}

SPRNode::SPRNode(std::string n)
{
    init(NULL, NULL, NULL, n, 0);
}

SPRNode::SPRNode(std::string n, int d)
{
    init(NULL, NULL, NULL, n, d);
}

SPRNode::SPRNode(SPRNode *lc, SPRNode *rc, SPRNode *p, std::string n, int d)
{
    init(lc, rc, p, n, d);
}

void SPRNode::init(SPRNode *lc, SPRNode *rc, SPRNode *p, std::string n, int d)
{
    //		this->lc = lc;
    //		this->rc = rc;
    this->p = p;
    this->name = std::string(n);
    this->twin = NULL;
    this->SPR_depth = d;
    this->pre_num = -1;
    this->edge_pre_start = -1;
    this->edge_pre_end = -1;
    this->component_number = -2;
    this->active_descendants = std::list <SPRNode *>();
    this->root_lcas = std::list <SPRNode *>();
    this->removable_descendants = std::list< std::list<SPRNode *>::iterator>();
    this->sibling_pair_loc = std::list<SPRNode *>::iterator();
    this->sibling_pair_status = 0;
    this->num_clustered_children = 0;
    this->forest = NULL;
    
    this->contracted_lc = NULL;
    this->contracted_rc = NULL;
    this->contracted = false;
    this->edge_protected = false;
    this->allow_sibling = true;
    this->lost_children = 0;
    this->max_merge_SPR_depth = -1;
    this->support = -1;
    this->support_normalization = -1;
    this->children = std::list<SPRNode *>();
    if (lc != NULL)
        add_child(lc);
    if (rc != NULL)
        add_child(rc);
}

SPRNode::SPRNode(const SPRNode &n)
{
    p = NULL;
    name = n.name.c_str();
    twin = n.twin;
    SPR_depth = n.SPR_depth;
    //		SPR_depth = 0;
    pre_num = n.pre_num;
    edge_pre_start = n.edge_pre_start;
    edge_pre_end = n.edge_pre_end;
    component_number = n.component_number;
    this->active_descendants = std::list <SPRNode *>();
    this->removable_descendants = std::list< std::list<SPRNode *>::iterator>();
    this->root_lcas = std::list <SPRNode *>();
    //sibling_pair_loc = n.sibling_pair_loc;
    //sibling_pair_status = n.sibling_pair_status;
    this->sibling_pair_loc = std::list<SPRNode *>::iterator();
    this->sibling_pair_status = 0;
    this->num_clustered_children = 0;
    this->forest = NULL;
    std::list<SPRNode *>::const_iterator c;
    this->children = std::list<SPRNode *>();
    for(c = n.children.begin(); c != n.children.end(); c++) {
        add_child(new SPRNode(**c));
    }
    
#ifdef COPY_CONTRACTED
    if (n.contracted_lc == NULL)
        contracted_lc = NULL;
    else
        contracted_lc = new SPRNode(*(n.contracted_lc), this);
    if (n.contracted_rc == NULL)
        contracted_rc = NULL;
    else
        contracted_rc = new SPRNode(*(n.contracted_rc), this);
    this->contracted = n.contracted;
#else
    this->contracted_lc = n.contracted_lc;
    this->contracted_rc = n.contracted_rc;
    this->contracted = n.contracted;
#endif
    this->edge_protected = n.edge_protected;
    this->allow_sibling = n.allow_sibling;
    this->lost_children = n.lost_children;
    this->max_merge_SPR_depth = n.max_merge_SPR_depth;
    this->support = n.support;
    this->support_normalization = n.support_normalization;
}

SPRNode::SPRNode(const SPRNode &n, SPRNode *parent)
{
    p = parent;
    name = n.name.c_str();
    twin = n.twin;
    if (p != NULL)
        SPR_depth = p->SPR_depth+1;
    else
        SPR_depth = n.SPR_depth;
    pre_num = n.pre_num;
    edge_pre_start = n.edge_pre_start;
    edge_pre_end = n.edge_pre_end;
    component_number = n.component_number;
    this->active_descendants = std::list <SPRNode *>();
    this->removable_descendants = std::list< std::list<SPRNode *>::iterator>();
    this->root_lcas = std::list <SPRNode *>();
    //sibling_pair_loc = n.sibling_pair_loc;
    //sibling_pair_status = n.sibling_pair_status;
    this->sibling_pair_loc = std::list<SPRNode *>::iterator();
    this->sibling_pair_status = 0;
    this->num_clustered_children = 0;
    this->forest = NULL;
    this->children = std::list<SPRNode *>();
    std::list<SPRNode *>::const_iterator c;
    for(c = n.children.begin(); c != n.children.end(); c++) {
        add_child(new SPRNode(**c));
    }
    
#ifdef COPY_CONTRACTED
    if (n.contracted_lc == NULL)
        contracted_lc = NULL;
    else
        contracted_lc = new SPRNode(*(n.contracted_lc), this);
    if (n.contracted_rc == NULL)
        contracted_rc = NULL;
    else
        contracted_rc = new SPRNode(*(n.contracted_rc), this);
    this->contracted = n.contracted;
#else
    this->contracted_lc = n.contracted_lc;
    this->contracted_rc = n.contracted_rc;
    this->contracted = n.contracted;
#endif
    this->edge_protected = n.edge_protected;
    this->allow_sibling = n.allow_sibling;
    this->lost_children = n.lost_children;
    this->max_merge_SPR_depth = n.max_merge_SPR_depth;
    this->support = n.support;
    this->support_normalization = n.support_normalization;
}

SPRNode::~SPRNode()
{
    std::list<SPRNode *>::iterator c = children.begin();
    while(c!= children.end()) {
        SPRNode *n = *c;
        c++;
        n->cut_parent();
    }
    cut_parent();
    active_descendants.clear();
    root_lcas.clear();
    removable_descendants.clear();
#ifdef COPY_CONTRACTED
    if (contracted_lc != NULL) {
        contracted_lc->delete_tree();
    }
    contracted_lc = NULL;
    if (contracted_rc != NULL) {
        contracted_rc->delete_tree();
    }
    contracted_rc = NULL;
#endif
}

// delete a subtree
void SPRNode::delete_tree()
{
    std::list<SPRNode *>::iterator c = children.begin();
    while(c!= children.end()) {
        SPRNode *n = *c;
        c++;
        n->delete_tree();
    }
#ifdef COPY_CONTRACTED
    if (contracted_lc != NULL) {
        contracted_lc->delete_tree();
    }
    contracted_lc = NULL;
    if (contracted_rc != NULL) {
        contracted_rc->delete_tree();
    }
    contracted_rc = NULL;
#endif
    delete this;
}

// cut edge between parent and child
// should really be cut_child, no deleting occurs
// TODO: is this useful? The only reason for this would be to
// be sure it works when the child is not correctly set
void SPRNode::delete_child(SPRNode *n)
{
    n->cut_parent();
}

// TODO: make sure this doesn't break things with >2 children
// add a child
void SPRNode::add_child(SPRNode *n)
{
    if (n->p != NULL)
        n->cut_parent();
    n->p_link = children.insert(children.end(),n);
    n->p = this;
    n->SPR_depth = SPR_depth+1;
    n->contracted = false;
}

// TODO: make sure this doesn't break things with >2 children
// add a child
void SPRNode::add_child_keep_SPR_depth(SPRNode *n)
{
    if (n->p != NULL)
        n->cut_parent();
    n->p_link = children.insert(children.end(),n);
    n->p = this;
    n->contracted = false;
}

// insert a child before the given sibling
void SPRNode::insert_child(SPRNode *sibling, SPRNode *n)
{
    if (n->p != NULL)
        n->cut_parent();
    n->p_link = children.insert(sibling->p_link, n);
    n->SPR_depth = SPR_depth+1;
    n->p = this;
    n->contracted = false;
}

// insert a child before the given sibling
void SPRNode::insert_child_keep_SPR_depth(SPRNode *sibling, SPRNode *n)
{
    if (n->p != NULL)
        n->cut_parent();
    n->p_link = children.insert(sibling->p_link, n);
    n->p = this;
    n->contracted = false;
}

SPRNode *SPRNode::set_twin(SPRNode *n)
{
    twin = n;
    return twin;
}

void SPRNode::set_name(std::string n)
{
    name = std::string(n);
}

int SPRNode::set_SPR_depth(int d)
{
    SPR_depth = d;
    return SPR_depth;
}

void SPRNode::fix_SPR_depths()
{
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++)
    {
        (*c)->SPR_depth = SPR_depth+1;
        (*c)->fix_SPR_depths();
    }
}

int SPRNode::set_preorder_number(int p)
{
    pre_num = p;
    return pre_num;
}

int SPRNode::set_edge_pre_start(int p)
{
    edge_pre_start= p;
    return edge_pre_start;
}

int SPRNode::set_edge_pre_end(int p)
{
    edge_pre_end= p;
    return edge_pre_end;
}

void SPRNode::copy_edge_pre_interval(SPRNode *n)
{
    if (n->edge_pre_start > -1)
    {
        edge_pre_start = n->edge_pre_start;
    }
    if (n->edge_pre_end > -1)
    {
        edge_pre_end = n->edge_pre_end;
    }
}

int SPRNode::set_component_number(int c)
{
    component_number = c;
    return component_number;
}

std::list<SPRNode *>& SPRNode::get_children()
{
    return children;
}

SPRNode *SPRNode::get_contracted_lc()
{
    return contracted_lc;
}

SPRNode *SPRNode::get_contracted_rc()
{
    return contracted_rc;
}

SPRNode *SPRNode::set_contracted_lc(SPRNode *n)
{
    contracted_lc = n;
    return contracted_lc;
}

SPRNode *SPRNode::set_contracted_rc(SPRNode *n)
{
    contracted_rc = n;
    return contracted_rc;
}

bool SPRNode::is_protected()
{
    return edge_protected;
}

bool SPRNode::is_contracted()
{
    return contracted;
}

void SPRNode::protect_edge()
{
    edge_protected = true;
}

void SPRNode::unprotect_edge()
{
    edge_protected = false;
}

void SPRNode::unprotect_subtree()
{
    unprotect_edge();
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++)
    {
        (*c)->unprotect_edge();
    }
}

void SPRNode::protect_supported_edges()
{
    if (support > 0)
        edge_protected = true;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++)
    {
        (*c)->protect_supported_edges();
    }
}

bool SPRNode::can_be_sibling()
{
    return allow_sibling;
}

void SPRNode::disallow_siblings()
{
    allow_sibling = false;
}

void SPRNode::disallow_siblings_subtree()
{
    allow_sibling = false;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++)
    {
        (*c)->disallow_siblings_subtree();
    }
}

void SPRNode::allow_siblings()
{
    allow_sibling = true;
}

void SPRNode::allow_siblings_subtree()
{
    allow_sibling = true;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++)
    {
        (*c)->allow_siblings_subtree();
    }
}

int SPRNode::num_lost_children()
{
    return lost_children;
}

int SPRNode::get_max_merge_SPR_depth()
{
    return max_merge_SPR_depth;
}

void SPRNode::set_max_merge_SPR_depth(int d)
{
    max_merge_SPR_depth = d;
}

int SPRNode::count_lost_children_subtree()
{
    int lost_children_count = lost_children;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++)
    {
        lost_children_count += (*c)->count_lost_children_subtree();
    }
    return lost_children_count;
}

int SPRNode::count_lost_subtree()
{
    int lost_children_count = (lost_children > 0) ? 1 : 0;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++)
    {
        lost_children_count += (*c)->count_lost_subtree();
    }
    return lost_children_count;
}

void SPRNode::lost_child()
{
    lost_children++;
}

void SPRNode::no_lost_children()
{
    lost_children = 0;
}

void SPRNode::no_lost_children_subtree()
{
    lost_children = 0;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->no_lost_children_subtree();
    }
}

double SPRNode::get_support()
{
    return support;
}

void SPRNode::set_support(double s)
{
    support = s;
}

void SPRNode::a_inc_support()
{
#pragma omp atomic
    support += 1;
}

void SPRNode::a_dec_support()
{
#pragma omp atomic
    support -= 1;
}

double SPRNode::get_support_normalization()
{
    return support_normalization;
}

void SPRNode::set_support_normalization(double s)
{
    support_normalization = s;
}

void SPRNode::a_inc_support_normalization()
{
#pragma omp atomic
    support_normalization += 1;
}

void SPRNode::a_dec_support_normalization()
{
#pragma omp atomic
    support_normalization -= 1;
}

void SPRNode::normalize_support()
{
    if (support_normalization != 0)
        support /= support_normalization;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++)
    {
        (*c)->normalize_support();
    }
}

int SPRNode::get_component_number()
{
    return component_number;
}

void SPRNode::increase_clustered_children()
{
    num_clustered_children++;
}

void SPRNode::decrease_clustered_children()
{
    num_clustered_children--;
}

void SPRNode::set_num_clustered_children(int c)
{
    num_clustered_children = c;
}

int SPRNode::get_num_clustered_children()
{
    return num_clustered_children;
}

std::list <SPRNode *> *SPRNode::get_active_descendants()
{
    return &active_descendants;
}

std::list <SPRNode *> *SPRNode::get_root_lcas()
{
    return &root_lcas;
}

int SPRNode::get_sibling_pair_status()
{
    return sibling_pair_status;
}

void SPRNode::set_sibling_pair_status(int s)
{
    sibling_pair_status = s;
}

void SPRNode::set_forest(Forest *f)
{
    forest = f;
}

void SPRNode::set_forest_rec(Forest *f)
{
    forest = f;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->set_forest_rec(f);
    }
}

Forest *SPRNode::get_forest()
{
    return forest;
}

std::list<std::list <SPRNode *>::iterator> *SPRNode::get_removable_descendants()
{
    return &removable_descendants;
}

void SPRNode::initialize_component_number(int value)
{
    component_number = value;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++)
    {
        (*c)->initialize_component_number(value);
    }
}

void SPRNode::initialize_active_descendants(std::list <SPRNode *> value)
{
    active_descendants = value;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++)
    {
        (*c)->initialize_active_descendants(value);
    }
}

void SPRNode::initialize_root_lcas(std::list <SPRNode *> value)
{
    root_lcas = value;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->initialize_root_lcas(value);
    }
}

void SPRNode::initialize_removable_descendants(std::list<std::list <SPRNode *>::iterator> value)
{
    removable_descendants = value;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++)
    {
        (*c)->initialize_removable_descendants(value);
    }
}


/* contract:
 * if this SPRNode has a parent and one child then contract it out
 * if this SPRNode has a parent and no child then contract it out
 *   parent will have one child so contract it as well.
 * if this SPRNode has no parent and one child then take the
 *	child's children and contract it out
 * return the first degree two parent found or NULL if there
 * was no contraction
 */
// TODO: check handling twins for interior SPRNodes
SPRNode *SPRNode::contract(bool remove)
{
    SPRNode *parent = p;
    SPRNode *child;
    SPRNode *ret = NULL;
    // contract out this SPRNode and give child to parent
    if (parent != NULL) {
        //cout << p->str_subtree() << endl;
        if (children.size() == 1) {
            child = children.front();
            if (this == parent->children.back()) {
                parent->add_child_keep_SPR_depth(child);
                child->set_SPR_depth(SPR_depth);
            }
            else {
                std::list<SPRNode *>::iterator sib = p_link;
                sib++;
                SPRNode *sibling = *sib;
                parent->insert_child_keep_SPR_depth(sibling, child);
                child->set_SPR_depth(SPR_depth);
            }
            child->copy_edge_pre_interval(this);
            if (edge_protected && !child->is_protected())
                child->protect_edge();
            cut_parent();
            if (remove)
                delete this;
            ret = parent;
        }
        else if (children.empty()) {
            cut_parent();
            ret = parent->contract(remove);
            if (remove)
                delete this;
            //this->fake_delete();
        }
        else
            ret = this;
        
    }
    // if no parent then take children of single child and remove it
    else {
        
        // dead component or singleton, will be cleaned up by the forest
        if (children.empty()) {
            //std::cout << get_name() << std::endl;
            if (str() == "")
                name = DEAD_COMPONENT;
        }
        if (children.size() == 1) {
            child = children.front();
            child->cut_parent();
            
            /* cluster hack - if we delete a cluster SPRNode then
             * we may try to use it later. This only happens once
             * per cluster so we can spend linear time to update
             * the forest
             */
            if (child->num_clustered_children > 0) {
                delete_child(child);
                if (remove)
                    delete this;
                //this->fake_delete();
                ret = child;
            }
            else {
                // if child is a leaf then get rid of this so we don't lose refs
                // problem: if the child is not c, then we want to copy
                // otherwise we don't
                // copy other parameters and join the twin
                //to this if the child is a label
                //SPRNode *new_lc = child->lchild();
                //SPRNode *new_rc = child->rchild();
                if (child->is_leaf()) {
                    if (child->get_twin() != NULL) {
                        set_twin(child->get_twin());
                        child->get_twin()->set_twin(this);
                    }
                    name = child->get_name().c_str();
                    //						name = child->str();
                }
                child->cut_parent();
                std::list<SPRNode *>::iterator c = child->children.begin();
                while(c!= child->children.end()) {
                    SPRNode *new_child = *c;
                    c++;
                    new_child->cut_parent();
                    add_child(new_child);
                }
                if (child->contracted_lc != NULL)
                    contracted_lc = child->contracted_lc;
                if (child->contracted_rc != NULL)
                    contracted_rc = child->contracted_rc;
                pre_num = child->get_preorder_number();
                if (remove) {
                    child->contracted_lc = NULL;
                    child->contracted_rc = NULL;
                    delete child;
                }
                ret = this;
            }
        }
    }
    
    return ret;
}

SPRNode *SPRNode::contract()
{
    return contract(false);
}

bool SPRNode::contract_sibling_pair()
{
    if (lchild() != NULL && lchild()->is_leaf()
        && rchild() != NULL && rchild()->is_leaf()) {
        std::string new_name = "(" + lchild()->str() + "," + rchild()->str() + ")";
        set_name(new_name);
        lchild()->cut_parent();
        rchild()->cut_parent();
        return true;
    }
    return false;
}

// TODO: binary only
bool SPRNode::contract_sibling_pair_undoable()
{
    if (lchild() != NULL && lchild()->is_leaf()
        && rchild() != NULL && rchild()->is_leaf())
    {
        SPRNode *lc = lchild();
        SPRNode *rc = rchild();
        contracted_lc = lc;
        contracted_rc = rc;
        rc->cut_parent();
        lc->cut_parent();
        contracted_lc->contracted = true;
        contracted_rc->contracted = true;
        edge_protected = false;
        return true;
    }
    return false;
}

/* contract_sibling_pair_undoable
 * works with multifurcating trees
 * returns NULL for no contract, otherwise returns the contracted
 * parent of the SPRNodes
 */
SPRNode *SPRNode::contract_sibling_pair_undoable(SPRNode *child1, SPRNode *child2)
{
    if (child1->parent() != this ||
        child2->parent() != this)
        return NULL;
    if (children.size() == 2) {
        contract_sibling_pair_undoable();
        return this;
    }
    else {
        SPRNode *new_child = new SPRNode();
        // buggy
        //			new_child->set_preorder_number(pre_num);
        if (child1->get_preorder_number() < child2->get_preorder_number()) {
            new_child->set_preorder_number(child1->get_preorder_number());
        }
        else {
            new_child->set_preorder_number(child2->get_preorder_number());
        }
        add_child(new_child);
        new_child->add_child(child1);
        new_child->add_child(child2);
        edge_protected = false;
        //new_child->contract_sibling_pair_undoable();
        return new_child;
    }
}

// TODO: binary only
void SPRNode::undo_contract_sibling_pair()
{
    // hacky, might hide problems
    if (contracted_lc != NULL)
        add_child(contracted_lc);
    if (contracted_rc != NULL)
        add_child(contracted_rc);
    contracted_lc = NULL;
    contracted_rc = NULL;
}

void SPRNode::fix_contracted_order()
{
    if (twin != NULL && twin->contracted_lc->twin != contracted_lc) {
        SPRNode *swap = contracted_lc;
        contracted_lc = contracted_rc;
        contracted_rc = swap;
    }
}

// cut the edge between this SPRNode and its parent
void SPRNode::cut_parent()
{
    if (p != NULL) {
        // TODO hacky: fix this to use a multi list for contractions
        if (!contracted) {
            p->children.erase(p_link);
        }
        else {
            if (p->contracted_lc == this)
                p->contracted_lc = NULL;
            if (p->contracted_rc == this)
                p->contracted_rc = NULL;
        }
        p = NULL;
        p_link = children.end();
    }
}

// caution: destructive
void SPRNode::contract_SPRNode()
{
    if (p == NULL || is_leaf())
        return;
    std::list<SPRNode *>::iterator c = children.begin();
    while(c!= children.end()) {
        SPRNode *n = *c;
        c++;
        p->add_child(n);
    }
    delete this;
}

SPRNode *SPRNode::parent()
{
    return p;
}

inline SPRNode *SPRNode::lchild()
{
    if (children.empty())
        return NULL;
    else
        return children.front();
}

SPRNode *SPRNode::rchild()
{
    if (children.empty())
        return NULL;
    std::list<SPRNode *>::iterator c = ++(children.begin());
    if (c == children.end())
        return NULL;
    else
        return *c;
}

SPRNode *SPRNode::get_twin()
{
    return twin;
}

int SPRNode::get_SPR_depth()
{
    return SPR_depth;
}

int SPRNode::get_preorder_number()
{
    return pre_num;
}

int SPRNode::get_edge_pre_start()
{
    return edge_pre_start;
}

int SPRNode::get_edge_pre_end()
{
    return edge_pre_end;
}

std::string SPRNode::str()
{
    std::string s = "";
    str_hlpr(&s);
    return s;
}

std::string SPRNode::get_name()
{
    return name;
}

void SPRNode::str_hlpr(std::string *s)
{
    if (!name.empty())
        *s += name.c_str();
    if (contracted_lc != NULL || contracted_rc != NULL) {
    *s += "(";
        if (contracted_lc != NULL) {
            contracted_lc->str_c_subtree_hlpr(s);
        }
        *s += ",";
        if (contracted_rc != NULL) {
            contracted_rc->str_c_subtree_hlpr(s);
        }
        *s += ")";
    }
}

std::string SPRNode::str_subtree()
{
    std::string s = "";
    str_subtree_hlpr(&s);
    return s;
}

void SPRNode::str_subtree_hlpr(std::string *s)
{
    str_hlpr(s);
    if (!is_leaf()) {
        *s += "(";
        std::list<SPRNode *>::iterator c;
        for(c = children.begin(); c != children.end(); c++) {
            if (c != children.begin())
                *s += ",";
            (*c)->str_subtree_hlpr(s);
            if ((*c)->parent() != this)
                std::cout << "#";
        }
        *s += ")";
    }
}

std::string SPRNode::str_support_subtree(bool allow_negative)
{
    std::string s = "";
    str_support_subtree_hlpr(&s, allow_negative);
    return s;
}

std::string SPRNode::str_support_subtree()
{
    return str_support_subtree(false);
}

void SPRNode::str_support_subtree_hlpr(std::string *s, bool allow_negative)
{
    str_hlpr(s);
    if (!is_leaf()) {
        *s += "(";
        std::list<SPRNode *>::iterator c;
        for(c = children.begin(); c != children.end(); c++) {
            if (c != children.begin())
                *s += ",";
            (*c)->str_support_subtree_hlpr(s, allow_negative);
            if ((*c)->parent() != this)
                std::cout << "#";
        }
        *s += ")";
        if (get_support() > -1 || allow_negative) {
            std::stringstream ss;
            ss << std::setprecision (2) << get_support();
            *s+= ss.str();
            if (get_support_normalization() > -1 || allow_negative) {
                std::stringstream ss;
                ss << "#";
                ss << get_support_normalization();
                *s+= ss.str();
            }
        }
    }
}

std::string SPRNode::str_edge_pre_interval_subtree()
{
    std::string s = "";
    str_edge_pre_interval_subtree_hlpr(&s);
    return s;
}

void SPRNode::str_edge_pre_interval_subtree_hlpr(std::string *s)
{
    str_hlpr(s);
    if (!is_leaf()) {
        *s += "(";
        std::list<SPRNode *>::iterator c;
        for(c = children.begin(); c != children.end(); c++) {
            if (c != children.begin())
                *s += ",";
            (*c)->str_edge_pre_interval_subtree_hlpr(s);
            if ((*c)->parent() != this)
                std::cout << "#";
        }
        *s += ")";
    }
    if (get_preorder_number() > -1) {
        std::stringstream ss;
        ss << ":";
        if (get_edge_pre_start() > -1) {
            ss << get_edge_pre_start();
            ss << ";";
        }
        ss << get_preorder_number();
        if (get_edge_pre_end() > -1) {
            ss << ";";
            ss << get_edge_pre_end();
        }
        *s+= ss.str();
    }
}

void SPRNode::str_c_subtree_hlpr(std::string *s)
{
    str_hlpr(s);
    if (!is_leaf()) {
        *s += "(";
        /*int cnt = 0;
        //std::cout << children.size() << std::endl;
        for(auto & c : children)
        {
            if(cnt) *s += ",";
            if(c != NULL)
            {
                //std::cout << c << std::endl;
                c->str_c_subtree_hlpr(s);
            }
        }*/
        std::list<SPRNode *>::iterator c;
        for(c = children.begin(); c != children.end(); c++) {
            std::cout << children.size() << std::endl;
            
                //std::cout << (*c)->get_name() << std::endl;
            if (c != children.begin())
                *s += ",";
            (*c)->str_c_subtree_hlpr(s);
            if ((*c)->parent() != this)
                std::cout << "#";
        }
        *s += ")";
    }
}

std::string SPRNode::str_subtree_twin()
{
    std::string s = "";
    str_subtree_twin_hlpr(&s);
    return s;
}

void SPRNode::str_subtree_twin_hlpr(std::string *s)
{
    *s += name.c_str();//str_hlpr(s);
    if (twin != NULL) {
        *s += "{";
        twin->str_subtree_hlpr(s);
        *s += "}";
    }
    if (!is_leaf()) {
        *s += "(";
        std::list<SPRNode *>::iterator c;
        for(c = children.begin(); c != children.end(); c++) {
            if (c != children.begin())
                *s += ",";
            (*c)->str_subtree_hlpr(s);
            if ((*c)->parent() != this)
                std::cout << "#";
        }
        *s += ")";
    }
}

void SPRNode::print()
{
    std::cout << name;
}

void SPRNode::print_subtree()
{
    std::cout << str_subtree();
    std::cout << std::endl;
}

void SPRNode::print_subtree_hlpr()
{
    std::cout << str_subtree();
}

void SPRNode::print_subtree_twin_hlpr()
{
    std::cout << str_subtree_twin();
}

bool SPRNode::is_leaf()
{
    return children.empty();
}

// TODO: binary only
bool SPRNode::is_sibling_pair()
{
    return (lchild() != NULL && lchild()->is_leaf()
            && rchild() != NULL && rchild()->is_leaf());
}

bool SPRNode::is_singleton()
{
    return (parent() == NULL && is_leaf());
}

// TODO: binary only
void SPRNode::find_sibling_pairs_hlpr(std::list<SPRNode *> *sibling_pairs)
{
    SPRNode *lchild = this->lchild();
    SPRNode *rchild = this->rchild();
    bool lchild_leaf = false;
    bool rchild_leaf = false;
    if (lchild != NULL) {
        if (lchild->is_leaf())
            lchild_leaf = true;
        else
            lchild->find_sibling_pairs_hlpr(sibling_pairs);
    }
    if (rchild != NULL) {
        if (rchild->is_leaf())
            rchild_leaf = true;
        else
            rchild->find_sibling_pairs_hlpr(sibling_pairs);
    }
    if (lchild_leaf && rchild_leaf) {
        sibling_pairs->push_back(lchild);
        sibling_pairs->push_back(rchild);
        //lchild->add_to_sibling_pairs(sibling_pairs, 1);
        //rchild->add_to_sibling_pairs(sibling_pairs, 2);
    }
}

// find the sibling pairs in this SPRNode's subtree
void SPRNode::append_sibling_pairs(std::list<SPRNode *> *sibling_pairs)
{
    find_sibling_pairs_hlpr(sibling_pairs);
}

// find the sibling pairs in this SPRNode's subtree
std::list<SPRNode *> *SPRNode::find_sibling_pairs()
{
    std::list<SPRNode *> *sibling_pairs = new std::list<SPRNode *>();
    find_sibling_pairs_hlpr(sibling_pairs);
    return sibling_pairs;
}

void SPRNode::find_leaves_hlpr(std::vector<SPRNode *> &leaves)
{
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        if ((*c)->is_leaf())
            leaves.push_back(*c);
        else
            (*c)->find_leaves_hlpr(leaves);
    }
    
}

// find the leaves in this SPRNode's subtree
std::vector<SPRNode *> SPRNode::find_leaves()
{
    std::vector<SPRNode *> leaves = std::vector<SPRNode *>();
    if (is_leaf())
        leaves.push_back(this);
    else
        find_leaves_hlpr(leaves);
    return leaves;
}

void SPRNode::find_interior_hlpr(std::vector<SPRNode *> &interior)
{
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        if (!(*c)->is_leaf()) {
            interior.push_back(*c);
            (*c)->find_interior_hlpr(interior);
        }
    }
    
}

// find the interior SPRNodes in this SPRNode's subtree
// does not include this SPRNode
std::vector<SPRNode *> SPRNode::find_interior()
{
    std::vector<SPRNode *> interior = std::vector<SPRNode *>();
    if (!is_leaf())
        find_interior_hlpr(interior);
    return interior;
}

void SPRNode::find_descendants_hlpr(std::vector<SPRNode *> &descendants)
{
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        descendants.push_back(*c);
        if (!(*c)->is_leaf()) {
            (*c)->find_descendants_hlpr(descendants);
        }
    }
    
}

// find the descendants SPRNodes in this SPRNode's subtree
// does not include this SPRNode
std::vector<SPRNode *> SPRNode::find_descendants()
{
    std::vector<SPRNode *> descendants = std::vector<SPRNode *>();
    if (!is_leaf())
        find_descendants_hlpr(descendants);
    return descendants;
}

bool SPRNode::contains_leaf(int number)
{
    if (stomini(get_name()) == number)
        return true;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        bool ans = (*c)->contains_leaf(number);
        if (ans == true)
            return true;
    }
    std::vector<SPRNode *> descendants = std::vector<SPRNode *>();
    return false;
}

// make twins point to this tree in this SPRNode's subtree
void SPRNode::resync()
{
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->resync();
    }
    if (twin != NULL)
        twin->set_twin(this);
}

// remove all twins
void SPRNode::unsync()
{
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->unsync();
    }
    twin = NULL;
}

// remove all interior twins
void SPRNode::unsync_interior()
{
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->unsync();
    }
    if(!is_leaf())
        twin = NULL;
}

void SPRNode::sync_af_twins()
{
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->sync_af_twins();
    }
    if (!is_leaf()) {
        twin = lchild()->get_twin()->parent();
        twin->set_twin(this);
    }
}

// find the root of this SPRNode's tree
SPRNode *SPRNode::find_root()
{
    SPRNode *root = this;
    while (root->parent() != NULL)
        root = root->parent();
    return root;
}

bool SPRNode::same_component(SPRNode *n, int &lca_SPR_depth, int &path_length)
{
    SPRNode *a = this;
    SPRNode *b = n;
    path_length = 0;
    while(a != b) {
        if ((b == NULL) || (a != NULL && a->get_SPR_depth() > b->get_SPR_depth()))
            a = a->parent();
        else
            b = b->parent();
        path_length++;
    }
    if (a == NULL) {
        path_length = -1;
        return false;
    }
    lca_SPR_depth = a->get_SPR_depth();
    return true;
}

bool SPRNode::same_component(SPRNode *n)
{
    int a,b;
    return same_component(n, a, b);
}

bool SPRNode::same_component(SPRNode *n, int &lca_SPR_depth)
{
    int a;
    return same_component(n, lca_SPR_depth, a);
}

void SPRNode::labels_to_numbers(std::map<std::string, int> *label_map, std::map<int, std::string> *reverse_label_map)
{
    if (name != "") {
        std::map<std::string, int>::iterator i = label_map->find(name);
        std::stringstream ss;
        if (i != label_map->end()) {
            ss << i->second;
            name = ss.str();
        }
        else {
            size_t num = label_map->size();
            ss << num;
            label_map->insert(make_pair(name, num));
            reverse_label_map->insert(make_pair(num, name));
            name = ss.str();
        }
    }
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->labels_to_numbers(label_map, reverse_label_map);
    }
    if (contracted_lc != NULL)
        contracted_lc->labels_to_numbers(label_map, reverse_label_map);
    if (contracted_rc != NULL)
        contracted_rc->labels_to_numbers(label_map, reverse_label_map);
}

void SPRNode::numbers_to_labels(std::map<int, std::string> *reverse_label_map)
{
    if (name != "") {
        std::string converted_name = "";
        std::string current_num = "";
        //string::iterator i = name.begin();
        size_t old_loc = 0;
        size_t loc = 0;
        while ((loc = name.find_first_of("0123456789", old_loc)) != std::string::npos) {
            converted_name.append(name.substr(old_loc, loc - old_loc));
            old_loc = loc;
            loc = name.find_first_not_of("0123456789", old_loc);
            std::string label = "";
            if (loc == std::string::npos)
                loc = name.size();
            label = name.substr(old_loc, loc - old_loc);
            std::map<int, std::string>::iterator j = reverse_label_map->find(atoi(label.c_str()));
            if (j != reverse_label_map->end()) {
                std::stringstream ss;
                ss << j->second;
                label = ss.str();
            }
            converted_name.append(label);
            old_loc = loc;
        }
        converted_name.append(name.substr(old_loc, name.size() - old_loc));
        name = converted_name;
    }
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->numbers_to_labels(reverse_label_map);
    }
    if (contracted_lc != NULL)
        contracted_lc->numbers_to_labels(reverse_label_map);
    if (contracted_rc != NULL)
        contracted_rc->numbers_to_labels(reverse_label_map);
}

void SPRNode::build_name_to_pre_map(std::map<std::string, int> *name_to_pre)
{
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->build_name_to_pre_map(name_to_pre);
    }
    if (is_leaf())
        name_to_pre->insert(make_pair(get_name(), get_preorder_number()));
}

void SPRNode::count_numbered_labels(std::vector<int> *label_counts)
{
    if (name != "") {
        int label = stomini(name);
        if (label_counts->size() <= label)
            label_counts->resize(label+1,0);
        (*label_counts)[label]++;
    }
    
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->count_numbered_labels(label_counts);
    }
    if (contracted_lc != NULL)
        contracted_lc->count_numbered_labels(label_counts);
    if (contracted_rc != NULL)
        contracted_rc->count_numbered_labels(label_counts);
}

void SPRNode::preorder_number()
{
    preorder_number(0);
}

int SPRNode::preorder_number(int next)
{
    set_preorder_number(next);
    next++;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        next = (*c)->preorder_number(next);
    }
    return next;
}

void SPRNode::edge_preorder_interval()
{
    edge_pre_start = pre_num;
    if (is_leaf()) {
        edge_pre_end = pre_num;
    }
    else {
        std::list<SPRNode *>::iterator c;
        for(c = children.begin(); c != children.end(); c++) {
            (*c)->edge_preorder_interval();
            if (edge_pre_end == -1 || (*c)->edge_pre_end > edge_pre_end)
                edge_pre_end = (*c)->edge_pre_end;
        }
    }
}

int SPRNode::size()
{
    int s  = 1;
    //std::cout << name << std::endl;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        //std::cout << (*c)->size() << std::endl;
        s += (*c)->size();
    }
    return s;
}

int SPRNode::size_using_prenum()
{
    if (is_leaf())
        return get_preorder_number();
    else
        return children.back()->size_using_prenum();
}

int SPRNode::max_SPR_depth()
{
    int d  = 0;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        int c_d = (*c)->max_SPR_depth();
        if (c_d > d)
            d = c_d;
    }
    return d+1;
}

size_t SPRNode::max_degree()
{
    size_t d = children.size();
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        size_t c_d = (*c)->max_degree();
        if (c_d > d)
            d = c_d;
    }
    return d;
}

// TODO: binary only
// these will potentially be removed

void SPRNode::add_to_front_sibling_pairs(std::list<SPRNode *> *sibling_pairs, int status)
{
    sibling_pairs->push_front(this);
    clear_sibling_pair(sibling_pairs);
    sibling_pair_status = status;
    sibling_pair_loc = sibling_pairs->begin();
}

void SPRNode::add_to_sibling_pairs(std::list<SPRNode *> *sibling_pairs, int status)
{
    sibling_pairs->push_back(this);
    clear_sibling_pair(sibling_pairs);
    sibling_pair_status = status;
    sibling_pair_loc = sibling_pairs->end();
    sibling_pair_loc--;
}

void SPRNode::remove_sibling_pair(std::list<SPRNode *> *sibling_pairs)
{
    if (sibling_pair_status > 0) {
        std::list<SPRNode *>::iterator loc = sibling_pair_loc;
        std::list<SPRNode *>::iterator sibling_loc = loc;
        if (sibling_pair_status == 1)
            sibling_loc++;
        else if (sibling_pair_status == 2)
            sibling_loc--;
        
        if (sibling_loc != sibling_pairs->end()) {
            SPRNode *old_sibling = *sibling_loc;
            old_sibling->sibling_pair_status = 0;
            sibling_pairs->erase(sibling_loc);
        }
        sibling_pairs->erase(loc);
        sibling_pair_status = 0;
    }
}

void SPRNode::clear_sibling_pair(std::list<SPRNode *> *sibling_pairs)
{
    if (sibling_pair_status > 0) {
        std::list<SPRNode *>::iterator loc = sibling_pair_loc;
        std::list<SPRNode *>::iterator sibling_loc = loc;
        if (sibling_pair_status == 1)
            sibling_loc++;
        else if (sibling_pair_status == 2)
            sibling_loc--;
        
        if (sibling_loc != sibling_pairs->end()) {
            SPRNode *old_sibling = *sibling_loc;
            old_sibling->sibling_pair_status = 0;
        }
        sibling_pair_status = 0;
    }
}


SPRNode *SPRNode::get_sibling()
{
    SPRNode *ret = NULL;
    if (p != NULL && p->children.size() > 1) {
        std::list<SPRNode *>::iterator s = p_link;
        if (p_link == p->children.begin())
            s++;
        else
            s--;
        ret = *s;
    }
    return ret;
}

SPRNode *SPRNode::get_right_sibling()
{
    SPRNode *ret = NULL;
    if (p != NULL && p->children.size() > 1) {
        std::list<SPRNode *>::iterator s = p_link;
        s++;
        if (s != children.end())
            ret = *s;
    }
    return ret;
}

// TODO: binary only
SPRNode *SPRNode::get_sibling(std::list<SPRNode *> *sibling_pairs)
{
    if (sibling_pair_status > 0) {
        std::list<SPRNode *>::iterator loc = sibling_pair_loc;
        std::list<SPRNode *>::iterator sibling_loc = loc;
        if (sibling_pair_status == 1)
            sibling_loc++;
        else if (sibling_pair_status == 2)
            sibling_loc--;
        return *sibling_loc;
    }
    else
        return NULL;
}

void SPRNode::set_sibling(SPRNode *sibling)
{
    if (sibling->sibling_pair_status > 0) {
        sibling_pair_loc = sibling->sibling_pair_loc;
        if (sibling->sibling_pair_status == 1)
            sibling_pair_loc++;
        else if (sibling->sibling_pair_status == 2)
            sibling_pair_loc--;
    }
}

void SPRNode::clear_sibling_pair_status()
{
    sibling_pair_status = 0;
}

void SPRNode::fix_parents()
{
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        if ((*c)->parent() != this) {
            (*c)->p = this;
            (*c)->p_link = c;
        }
        (*c)->fix_parents();
    }
}

void SPRNode::left_rotate()
{
    if (lchild()->lchild() != NULL) {
        SPRNode *new_lc = lchild()->lchild();
        SPRNode *new_rc = lchild();
        SPRNode *new_rc_lc = lchild()->rchild();
        SPRNode *new_rc_rc = rchild();
        
        new_lc->cut_parent();
        new_rc->cut_parent();
        new_rc_lc->cut_parent();
        new_rc_rc->cut_parent();
        add_child(new_lc);
        add_child(new_rc);
        rchild()->add_child(new_rc_lc);
        rchild()->add_child(new_rc_rc);
    }
}

void SPRNode::right_rotate()
{
    if (rchild()->lchild() != NULL) {
        SPRNode *new_lc = rchild()->lchild();
        SPRNode *new_rc = rchild();
        SPRNode *new_rc_lc = rchild()->rchild();
        SPRNode *new_rc_rc = lchild();
        
        new_lc->cut_parent();
        new_rc->cut_parent();
        new_rc_lc->cut_parent();
        new_rc_rc->cut_parent();
        add_child(new_lc);
        add_child(new_rc);
        rchild()->add_child(new_rc_lc);
        rchild()->add_child(new_rc_rc);
    }
}

void SPRNode::next_rooting()
{
    if (lchild()->lchild() != NULL)
        left_rotate();
    else if (rchild()->lchild() != NULL)
        right_rotate();
    else
        return;
    if (lchild()->pre_num < rchild()->pre_num && (lchild()->pre_num != 1 || rchild()->pre_num != 2 || !lchild()->is_leaf()))
        next_rooting();
}

/* reroot
 * rerooots the tree between new_lc and its parent. Maintains this
 * SPRNode as the root
 * assumes that new_lc is a descendant of this SPRNode
 * assumes this SPRNode is the root and has two children
 * other SPRNodes may be multifurcating
 */
void SPRNode::reroot(SPRNode *new_lc)
{
    SPRNode *new_rc = new_lc->parent();
    if (new_rc == this || new_rc == NULL) {
        return;
    }
    SPRNode *prev = new_rc;
    SPRNode *next = new_rc->parent();
    //	SPRNode *old_lc = lchild();
    //	SPRNode *old_rc_rc = rchild();
    new_lc->cut_parent();
    new_rc->cut_parent();
    while(next != NULL) {
        SPRNode *current = next;
        next = current->parent();
        prev->add_child(current);
        prev = current;
    }
    SPRNode *root = prev;
    while(root->get_children().size() > 0)
        root->parent()->add_child(root->lchild());
    //	root->parent()->add_child(root->lchild());
    root->cut_parent();
    root->add_child(new_lc);
    root->add_child(new_rc);
}

// make the root binay again
void SPRNode::fixroot()
{
    if (get_children().size() > 2) {
        SPRNode *new_lc = new SPRNode();
        while(get_children().size() > 1) {
            new_lc->add_child(get_children().back());
        }
        add_child(new_lc);
    }
}


SPRNode *SPRNode::expand_parent_edge(SPRNode *n)
{
    if (p != NULL) {
        SPRNode *old_p = p;
        cut_parent();
        SPRNode *new_p = new SPRNode();
        new_p->add_child(n);
        old_p->add_child(new_p);
        return p;
    }
    else {
        SPRNode *new_child = new SPRNode(name);
        SPRNode *old_lc = lchild();
        SPRNode *old_rc = rchild();
        old_lc->cut_parent();
        old_rc->cut_parent();
        
        new_child->add_child(old_lc);
        new_child->add_child(old_rc);
        new_child->contracted_lc = contracted_lc;
        if (contracted_lc != NULL)
            contracted_lc->p = new_child;
        new_child->contracted_rc = contracted_rc;
        if (contracted_rc != NULL)
            contracted_rc->p = new_child;
        name = "";
        contracted_lc = NULL;
        contracted_rc = NULL;
        add_child(new_child);
        return this;
    }
}

SPRNode *SPRNode::undo_expand_parent_edge() {
    if (p != NULL) {
        SPRNode *old_p = p;
        cut_parent();
        SPRNode *child = children.front();
        child->cut_parent();
        old_p->add_child(child);
        return this;
    }
    else {
        SPRNode *child = children.front();
        name = child->name;
        SPRNode *new_lc = child->lchild();
        SPRNode *new_rc = child->rchild();
        new_lc->cut_parent();
        new_rc->cut_parent();
        contracted_lc = child->contracted_lc;
        if (contracted_lc != NULL)
            contracted_lc->p = this;
        contracted_rc = child->contracted_rc;
        if (contracted_rc != NULL)
            contracted_rc->p = this;
        child->contracted_lc = NULL;
        child->contracted_rc = NULL;
        child->name = "";
        add_child(new_lc);
        add_child(new_rc);
        return child;
    }
}

/*  apply an SPR operation to move this subtree to be a
 *	sibling of new_sibling
 *
 *  Note: problems will occur if new_sibling is a descendant of this
 *  Note: does nothing if this is the root
 *
 *	Returns a SPRNode that will reverse the spr (old_sibling unless it
 *		was moved to maintain a root, in which case the root is returned,
 *		NULL if no spr occured)
 *
 * Note: binary only
 */
SPRNode *SPRNode::spr(SPRNode *new_sibling, int &which_child) {
    SPRNode *reverse;
    int prev_child_loc = 0;
    if (p == NULL || new_sibling == NULL)
        return NULL;
    SPRNode *old_sibling = get_sibling();
    if (old_sibling == new_sibling)
        return NULL;
    SPRNode *grandparent = p->p;
    if (p->lchild() == this)
        prev_child_loc = 1;
    else
        prev_child_loc = 2;
    // Prune
    if (grandparent != NULL) {
        bool leftc = false;
        if (old_sibling->parent() == grandparent->lchild())
            leftc = true;
        
        old_sibling->cut_parent();
        //p->delete_child(old_sibling);
        p->cut_parent();
        //		grandparent->delete_child(p);
        if (leftc && !grandparent->is_leaf()) {
            SPRNode *ns = grandparent->children.front();
            grandparent->insert_child(ns,old_sibling);
        }
        else
            grandparent->add_child(old_sibling);
        
        reverse = old_sibling;
    }
    else {
        if (old_sibling->is_leaf())
            return NULL;
        SPRNode *root = p;
        bool leftc = false;
        if (root->lchild() == this)
            leftc = true;
        SPRNode *lc = old_sibling->lchild();
        SPRNode *rc = old_sibling->rchild();
        root->delete_child(this);
        root->delete_child(old_sibling) ;
        if (lc != NULL)
            root->add_child(lc);
        if (rc != NULL)
            root->add_child(rc);
        
        //old_sibling->delete_child(old_sibling->lchild());
        //old_sibling->delete_child(old_sibling->rchild());
        if (leftc) {
            if (old_sibling->is_leaf())
                old_sibling->add_child(this);
            else {
                SPRNode *new_sibling = old_sibling->children.front();
                old_sibling->insert_child(new_sibling,this);
            }
        }
        else
            old_sibling->add_child(this);
        reverse = root;
    }
    
    
    // Regraft
    if (new_sibling->p != NULL) {
        grandparent = new_sibling->p;
        //		grandparent->delete_child(new_sibling);
        bool leftc = false;
        if (new_sibling == grandparent->lchild())
            leftc = true;
        if (leftc && !grandparent->is_leaf()) {
            SPRNode *new_sibling = grandparent->children.front();
            grandparent->insert_child(new_sibling,p);
        }
        else {
            grandparent->add_child(p);
        }
        
        if (which_child == 1)
            p->add_child(new_sibling);
        else {
            SPRNode *ns = p->children.front();
            p->insert_child(ns,new_sibling);
        }
        
    }
    else {
        SPRNode *root = new_sibling;
        new_sibling = p;
        // TODO: still broken
        SPRNode *lc = root->lchild();
        SPRNode *rc = root->rchild();
        p->add_child(lc);
        p->add_child(rc);
        //		p->delete_child(this);
        // problem here
        if (which_child == 0)
            which_child = prev_child_loc;
        if (which_child == 1) {
            root->add_child(this);
            root->add_child(new_sibling);
        }
        else {
            root->add_child(new_sibling);
            root->add_child(this);
        }
        
        
    }
    
    which_child = prev_child_loc;
    return reverse;
}

void SPRNode::spr(SPRNode *new_sibling)
{
    int na = 0;
    spr(new_sibling, na);
}

void SPRNode::find_descendant_counts_hlpr(std::vector<int> *dc)
{
    int num_descendants = 0;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->find_descendant_counts_hlpr(dc);
        num_descendants += (*dc)[(*c)->get_preorder_number()];
        num_descendants += 1;
    }
    if (dc->size() <= get_preorder_number())
        dc->resize(get_preorder_number() + 1, -1);
    (*dc)[get_preorder_number()] = num_descendants;
}

std::vector<int> *SPRNode::find_descendant_counts()
{
    std::vector<int> *dc = new std::vector<int>();
    find_descendant_counts_hlpr(dc);
    return dc;
}

void SPRNode::find_leaf_counts_hlpr(std::vector<int> *lc)
{
    int num_leaves = 0;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        (*c)->find_leaf_counts_hlpr(lc);
        num_leaves += (*lc)[(*c)->get_preorder_number()];
    }
    if (is_leaf())
        num_leaves++;
    if (lc->size() <= get_preorder_number())
        lc->resize(get_preorder_number() + 1, -1);
    (*lc)[get_preorder_number()] = num_leaves;
}

std::vector<int> *SPRNode::find_leaf_counts()
{
    std::vector<int> *lc = new std::vector<int>();
    find_leaf_counts_hlpr(lc);
    return lc;
}

SPRNode *SPRNode::find_median()
{
    std::vector<int> *dc = find_descendant_counts();
    return find_median_hlpr(dc, (*dc)[get_preorder_number()] / 2);
}

SPRNode *SPRNode::find_subtree_of_size(double percentage)
{
    std::vector<int> *dc = find_descendant_counts();
    return find_median_hlpr(dc, (int)((*dc)[get_preorder_number()] * percentage));
}

SPRNode *SPRNode::find_median_hlpr(std::vector<int> *dc, int target_size)
{
    SPRNode *largest_child_subtree = NULL;
    int lcs_size = 0;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        int cs_size = (*dc)[(*c)->get_preorder_number()] + 1;
        if (cs_size > lcs_size) {
            largest_child_subtree = *c;
            lcs_size = cs_size;
        }
    }
    if (lcs_size > target_size)
        return largest_child_subtree->find_median_hlpr(dc, target_size);
    else
        return this;
}

int SPRNode::any_leaf_preorder_number()
{
    if (is_leaf()) {
        if (contracted_lc != NULL)
            return contracted_lc->any_leaf_preorder_number();
        else return get_preorder_number();
    }
    else
        return (*children.begin())->any_leaf_preorder_number();
}

SPRNode *SPRNode::find_by_prenum(int prenum)
{
    if (prenum == get_preorder_number())
        return this;
    SPRNode *search_child = NULL;
    int best_prenum = -1;
    std::list<SPRNode *>::iterator c;
    for(c = children.begin(); c != children.end(); c++) {
        int p = (*c)->get_preorder_number();
        if (p == prenum)
            return *c;
        else if (p < prenum && p > best_prenum) {
            best_prenum = p;
            search_child = *c;
        }
    }
    if (search_child != NULL)
        return search_child->find_by_prenum(prenum);
    else
        return NULL;
}

// expand all contracted SPRNodes of a subtree starting at n
void SPRNode::expand_contracted_SPRNodes()
{
    if (is_leaf()) {
        if (contracted_lc != NULL) {
            add_child(contracted_lc);
            contracted_lc = NULL;
        }
        if (contracted_rc != NULL) {
            add_child(contracted_rc);
            contracted_rc = NULL;
        }
    }
    std::list<SPRNode *>::iterator c;
    for(c = get_children().begin(); c != get_children().end(); c++) {
        (*c)->expand_contracted_SPRNodes();
    }
}

int SPRNode::get_name_num() {
    return atoi(get_name().c_str());
}

// build a tree from a newick string
SPRNode *spr_building_tree(std::string s) {
    return spr_building_tree(s, 0, NULL);
}
SPRNode *spr_building_tree(std::string s, int start_SPR_depth) {
    return spr_building_tree(s, start_SPR_depth, NULL);
}
SPRNode *spr_building_tree(std::string s, std::set<std::string, StringCompare> *include_only) {
    return spr_building_tree(s, 0, include_only);
}
SPRNode *spr_building_tree(std::string s, int start_SPR_depth, std::set<std::string, StringCompare> *include_only) {
    if (s == "")
        return new SPRNode();
    SPRNode *dummy_head = new SPRNode("p", start_SPR_depth-1);
    bool valid = true;
    spr_building_tree_helper(0, s, dummy_head, valid, include_only);
    SPRNode *head = dummy_head->lchild();
    if (valid && head != NULL) {
        delete dummy_head;
        return head;
    }
    else {
        if (head != NULL)
            head->delete_tree();
        return dummy_head;
    }
    
}

// spr_building_tree recursive helper function
size_t spr_building_tree_helper(size_t start, const std::string& s, SPRNode *parent,
                                bool &valid, std::set<std::string, StringCompare> *include_only)
{
    size_t loc = s.find_first_of("(,)", start);
    if (loc == std::string::npos) {
        std::string name = s.substr(start, s.size() - start);
        size_t name_end = name.find(':');
        if (name_end != std::string::npos)
            name = name.substr(0, name_end);
        if (include_only == NULL ||
            include_only->find(name) != include_only->end()) {
            SPRNode *node = new SPRNode(name);
            parent->add_child(node);
        }
        loc = s.size()-1;
        return loc;
    }
    while(s[start] == ' ' || s[start] == '\t')
        start++;
    size_t end = loc;
    while(s[end] == ' ' || s[end] == '\t')
        end--;
    std::string name = s.substr(start, end - start);
    size_t name_end = name.find(':');
    if (name_end != std::string::npos)
        name = name.substr(0, name_end);
    SPRNode *node = NULL;
    if (include_only == NULL ||
        include_only->find(name) != include_only->end()) {
        node = new SPRNode(name);
        parent->add_child(node);
    }
    
    int count = 1;
    if (s[loc] == '(') {
        loc = spr_building_tree_helper(loc + 1, s, node, valid, include_only);
        while(s[loc] == ',') {
            loc = spr_building_tree_helper(loc + 1, s, node, valid, include_only);
            count++;
        }
        //			int loc_check = s.find_first_of("(,)", loc);
        //			if (loc_check != string::npos &&
        //					s[loc_check] == ','
        if (s[loc] != ')'
            || (false && count > 2)) {
            //|| (IGNORE_MULTI && count > 2)) {
            valid = false;
            return s.size()-1;
        }
        // TODO: get the support values here (and branch lengths?)
        // contract_SPRNode() if support is less than a threshold
        loc++;
        if (s[loc-1] == ')') {
            size_t numc = node->get_children().size();
            bool contracted = false;
            size_t next = s.find_first_of(",)", loc);
            if (next != std::string::npos) {
                if (next > loc && 0.0 > 0) {
                    //if (next > loc && REQUIRED_SUPPORT > 0) {
                    std::string info = s.substr(loc, next - loc);
                    if (info[0] != ':') {
                        double support = atof(info.c_str());
                        //							cout << "support=" << support << endl;
                        //if (support < REQUIRED_SUPPORT && numc > 0) {
                        if (support < 0.0 && numc > 0) {
                            node->contract_SPRNode();
                            contracted = true;
                        }
                    }
                }
                loc=next;
            }
            if (!contracted) {
                if (numc == 1)
                    node->contract_SPRNode();
                else if (numc == 0 && name == "") {
                    node->cut_parent();
                    delete node;
                }
            }
        }
    }
    return loc;
}

// swap two SPRNodes
void swap(SPRNode **a, SPRNode **b) {
    SPRNode *temp = *a;
    *a = *b;
    *b = temp;
}

int stomini(std::string s) {
    //	cout << "stomini" << endl;
    std::string number_characters = "+-0123456789";
    int min = INT_MAX;
    std::string current = "";
    for(int i = 0; i < s.size(); i++) {
        if (number_characters.find(s[i]) != std::string::npos) {
            current += s[i];
        }
        else if (current.size() > 0) {
            int num = atoi(current.c_str());
            if (num < min)
                min = num;
            current = "";
        }
    }
    if (current.size() > 0) {
        int num = atoi(current.c_str());
        if (num < min)
            min = num;
        current = "";
    }
    //	cout << "returning " << min << endl;
    return min;
}

// assumes that an unrooted tree is represented with a 3-way multifurcation
std::string SPR_root(std::string s) {
    //	cout << "root(string s)" << endl;
    //	cout << s << endl;
    std::string r = "";
    int SPR_depth = 0;
    int first_c = -1;
    int second_c = -1;
    int last_bracket = -1;
    for(int i = 0; i < s.size(); i++) {
        if (s[i] == '(')
            SPR_depth++;
        else if (s[i] == ')') {
            SPR_depth--;
            last_bracket = i;
        }
        else if (SPR_depth == 1 && s[i] == ',') {
            if (first_c == -1)
                first_c = i;
            else if  (second_c == -1)
                second_c = i;
        }
    }
    if (second_c == -1 || last_bracket == -1)
        return s;
    else {
        r.append(s.substr(0,first_c+1));
        r.append("(");
        r.append(s.substr(first_c+1,last_bracket-first_c));
        r.append(")");
        r.append(s.substr(last_bracket+1,std::string::npos));
    }
    //	cout << r << endl;
    return r;
}






































