#ifndef UndoMachine_hpp
#define UndoMachine_hpp

#include "SPRNode.hpp"
#include "Forest.hpp"
#include "SiblingPair.hpp"
#include <list>
#include <typeinfo>

class SPRNode;

class Undoable{
public:
    virtual ~Undoable() {}
    virtual void undo() = 0;
};

class UndoMachine{
public:
    std::list<Undoable *> events;
    int size;
    
public:
    UndoMachine();
    void add_event(Undoable *event);
    void insert_event(std::list<Undoable *>::iterator i, Undoable *event);
    std::list<Undoable *>::iterator get_bookmark();
    void undo();
    void undo(int num);
    void undo_all();
    void undo_to(int to);
    void clear_to(int to);
    int num_events();
};

class AddRho : public Undoable {
public:
    Forest *F;
    AddRho(Forest *f) {
        F = f;
    }
    void undo() {
        F->set_rho(false);
        F->get_component((int)F->num_components()-1)->delete_tree();
        F->erase_components((int)F->num_components()-1,(int)F->num_components());
    }
};

class AddComponent : public Undoable {
public:
    Forest *F;
    AddComponent(Forest *f) {
        F = f;
    }
    void undo() {
        F->erase_components((int)F->num_components()-1,(int)F->num_components());
    }
};

class AddComponentToFront : public Undoable {
public:
    Forest *F;
    AddComponentToFront(Forest *f) {
        F = f;
    }
    void undo() {
        F->erase_components(0,1);
    }
};


class CutParent : public Undoable {
public:
    SPRNode *child;
    SPRNode *parent;
    int branch;
    int SPR_depth;
    CutParent(SPRNode *c) {
        child = c;
        parent = c->parent();
        branch = 0;
        SPR_depth = c->get_SPR_depth();
        if (parent != NULL) {
            if (parent->lchild() == child)
                branch = 1;
            else
                branch = 2;
        }
    }
    
    void undo() {
        if (branch == 1) {
            if (parent->is_leaf())
                parent->add_child(child);
            else
                parent->insert_child(parent->get_children().front(), child);
        }
        else if (branch == 2)
            parent->add_child(child);
        child->set_SPR_depth(SPR_depth);
    }
};

class ClearSiblingPair : public Undoable {
public:
    SPRNode *a;
    SPRNode *c;
    int a_status;
    int c_status;
    ClearSiblingPair(SPRNode *x, SPRNode *y) {
        if (x->get_sibling_pair_status() == 1 ||
            y->get_sibling_pair_status() == 2) {
            a = x;
            c = y;
        }
        else {
            a = y;
            c = x;
        }
        a_status = a->get_sibling_pair_status();
        c_status = c->get_sibling_pair_status();
    }
    
    void undo() {
        a->set_sibling_pair_status(1);
        c->set_sibling_pair_status(2);
        if (a_status == 0) {
            c->set_sibling(a);
        }
        else if (c_status == 0) {
            a->set_sibling(c);
        }
        
    }
};

class PopClearedSiblingPair : public Undoable {
public:
    SPRNode *a;
    SPRNode *c;
    int a_status;
    int c_status;
    
    std::list<SPRNode *> *sibling_pairs;
    PopClearedSiblingPair(SPRNode *x, SPRNode *y, std::list<SPRNode *> *s) {
        sibling_pairs = s;
        a = x;
        c = y;
    }
    
    void undo() {
        sibling_pairs->push_back(c);
        sibling_pairs->push_back(a);
    }
};

class PopSiblingPair : public Undoable {
public:
    SPRNode *a;
    SPRNode *c;
    std::list<SPRNode *> *sibling_pairs;
    PopSiblingPair(SPRNode *x, SPRNode *y, std::list<SPRNode *> *s) {
        sibling_pairs = s;
        a = x;
        c = y;
    }
    
    void undo() {
        sibling_pairs->push_back(c);
        sibling_pairs->push_back(a);
    }
};

class ContractSiblingPair : public Undoable {
public:
    SPRNode *node;
    int c1_SPR_depth;
    int c2_SPR_depth;
    bool binary_node;
    bool node_protected;
    ContractSiblingPair(SPRNode *n) {
        init(n);
        if (n->is_protected())
            node_protected = true;
        else
            node_protected = false;
    }
    ContractSiblingPair(SPRNode *n, SPRNode *child1, SPRNode *child2,
                        UndoMachine *um) {
        if (n->get_children().size() == 2)
            init(n);
        else {
            um->add_event(new CutParent(child1));
            um->add_event(new CutParent(child2));
            binary_node = false;
            node = n;
        }
        if (n->is_protected())
            node_protected = true;
        else
            node_protected = false;
    }
    
    void init(SPRNode *n) {
        node = n;
        if (n->lchild() != NULL)
            c1_SPR_depth = n->lchild()->get_SPR_depth();
        else
            c1_SPR_depth = -1;
        if (n->rchild() != NULL)
            c2_SPR_depth = n->rchild()->get_SPR_depth();
        else
            c2_SPR_depth = -1;
        binary_node = true;
    }
    
    void undo() {
        if (binary_node) {
            node->undo_contract_sibling_pair();
            if (c1_SPR_depth > -1)
                node->lchild()->set_SPR_depth(c1_SPR_depth);
            if (c2_SPR_depth > -1)
                node->rchild()->set_SPR_depth(c2_SPR_depth);
        }
        if (node_protected)
            node->protect_edge();
    }
};

class AddToFrontSiblingPairs : public Undoable {
public:
    std::list<SPRNode *> *sibling_pairs;
    AddToFrontSiblingPairs(std::list<SPRNode *> *s) {
        sibling_pairs = s;
    }
    void undo() {
        if (!sibling_pairs->empty()) {
            sibling_pairs->pop_front();
            sibling_pairs->pop_front();
        }
    }
};

class AddToSiblingPairs : public Undoable {
public:
    std::list<SPRNode *> *sibling_pairs;
    AddToSiblingPairs(std::list<SPRNode *> *s) {
        sibling_pairs = s;
    }
    void undo() {
        if (!sibling_pairs->empty()) {
            sibling_pairs->pop_back();
            sibling_pairs->pop_back();
        }
    }
};

class AddToSetSiblingPairs : public Undoable {
public:
    std::set<SiblingPair> *sibling_pairs;
    SiblingPair pair;
    AddToSetSiblingPairs(std::set<SiblingPair> *sp, SiblingPair p) {
        sibling_pairs = sp;
        //pair = SiblingPair(p->a,p->c);
        pair = p;
    }
    void undo() {
        if (!sibling_pairs->empty()) {
            sibling_pairs->erase(pair);
        }
    }
};

class RemoveSetSiblingPairs : public Undoable {
public:
    std::set<SiblingPair> *sibling_pairs;
    SiblingPair pair;
    RemoveSetSiblingPairs(std::set<SiblingPair> *sp, SiblingPair p) {
        sibling_pairs = sp;
        pair = p;
    }
    void undo() {
        sibling_pairs->insert(pair);
    }
};

class AddInSiblingPairs : public Undoable {
public:
    std::list<SPRNode *> *sibling_pairs;
    int pos;
    AddInSiblingPairs(std::list<SPRNode *> *s, int p) {
        sibling_pairs = s;
        pos = p;
    }
    void undo() {
        if (!sibling_pairs->empty()) {
            std::list<SPRNode *>::iterator c = sibling_pairs->begin();
            for(int i = 0; i <= pos && c != sibling_pairs->end(); i++) {
                c++;
            }
            if (c != sibling_pairs->end()) {
                std::list<SPRNode *>::iterator rem = c;
                c++;
                sibling_pairs->erase(rem);
                rem = c;
                c++;
                sibling_pairs->erase(rem);
            }
        }
    }
};

class SetTwin : public Undoable {
public:
    SPRNode *node;
    SPRNode *twin;
    
    SetTwin(SPRNode *n) {
        node = n;
        twin = n->get_twin();
    }
    
    void undo() {
        node->set_twin(twin);
    }
};

class ChangeName : public Undoable {
public:
    SPRNode *node;
    std::string name;
    
    ChangeName(SPRNode *n) {
        node = n;
        name = n->get_name();//str();
        //			name = n->str();
    }
    
    void undo() {
        node->set_name(name);
    }
};

class ChangeEdgePreInterval : public Undoable {
public:
    SPRNode *node;
    int start;
    int end;
    
    ChangeEdgePreInterval(SPRNode *n) {
        start = n->get_edge_pre_start();
        end = n->get_edge_pre_end();
        node = n;
    }
    
    void undo() {
        node->set_edge_pre_start(start);
        node->set_edge_pre_end(end);
    }
};

class ChangePreNum : public Undoable {
public:
    SPRNode *node;
    int prenum;
    
    ChangePreNum(SPRNode *n) {
        prenum = n->get_preorder_number();
        node = n;
    }
    
    void undo() {
        node->set_preorder_number(prenum);
    }
};

class ChangeRightChild : public Undoable {
public:
    SPRNode *node;
    SPRNode *rchild;
    int rchild_SPR_depth;
    
    ChangeRightChild(SPRNode *n) {
        node = n;
        rchild = n->rchild();
        if (rchild != NULL)
            rchild_SPR_depth = rchild->get_SPR_depth();
    }
    
    void undo() {
        if (rchild != NULL) {
            node->add_child(rchild);
            rchild->set_SPR_depth(rchild_SPR_depth);
        }
        else
            if (node->rchild() != NULL)
                node->rchild()->cut_parent();
    }
};

class ChangeLeftChild : public Undoable {
public:
    SPRNode *node;
    SPRNode *lchild;
    int lchild_SPR_depth;
    
    ChangeLeftChild(SPRNode *n) {
        node = n;
        lchild = n->lchild();
        if (lchild != NULL)
            lchild_SPR_depth = lchild->get_SPR_depth();
    }
    
    void undo() {
        if (lchild != NULL) {
            //SPRNode->add_child_keep_SPR_depth(lchild);
            node->add_child(lchild);
            lchild->set_SPR_depth(lchild_SPR_depth);
        }
        else
            if (node->lchild() != NULL)
                node->lchild()->cut_parent();
    }
};

class AddChild : public Undoable {
public:
    SPRNode *child;
    int SPR_depth;
    
    AddChild(SPRNode *c) {
        child = c;
        if (c != NULL)
            SPR_depth = c->get_SPR_depth();
    }
    
    void undo() {
        if (child != NULL) {
            child->cut_parent();
            child->set_SPR_depth(SPR_depth);
        }
    }
};

class AddContractedLC : public Undoable {
public:
    SPRNode *node;
    
    AddContractedLC(SPRNode *n) {
        node = n;
    }
    
    void undo() {
        node->set_contracted_lc(NULL);
    }
};

class AddContractedRC : public Undoable {
public:
    SPRNode *node;
    
    AddContractedRC(SPRNode *n) {
        node = n;
    }
    
    void undo() {
        node->set_contracted_rc(NULL);
    }
};

class CreateSPRNode : public Undoable {
public:
    SPRNode *node;
    
    CreateSPRNode(SPRNode *n) {
        node = n;
    }
    
    void undo() {
        if (node != NULL)
            delete node;
    }
};

class ProtectEdge : public Undoable {
public:
    SPRNode *node;
    
    ProtectEdge(SPRNode *n) {
        node = n;
    }
    
    void undo() {
        if (node != NULL)
            node->unprotect_edge();
    }
};

class UnprotectEdge : public Undoable {
public:
    SPRNode *node;
    
    UnprotectEdge(SPRNode *n) {
        node = n;
    }
    
    void undo() {
        if (node != NULL)
            node->protect_edge();
    }
};

class ListPushBack : public Undoable {
public:
    std::list<SPRNode *> *l;
    
    ListPushBack(std::list<SPRNode *> *x) {
        l = x;
    }
    
    void undo() {
        l->pop_back();
    }
};

class ListPopBack : public Undoable {
public:
    std::list<SPRNode *> *l;
    SPRNode *node;
    
    ListPopBack(std::list<SPRNode *> *x) {
        l = x;
        if (!l->empty())
            node = l->back();
        else
            node = NULL;
    }
    
    void undo() {
        if (node != NULL)
            l->push_back(node);
    }
};

void ContractEvent(UndoMachine *um, SPRNode *n, std::list<Undoable *>::iterator
                   bookmark);

void ContractEvent(UndoMachine *um, SPRNode *n);
    
    
    




#endif /* UndoMachine_hpp */
