#include "UndoMachine.hpp"

UndoMachine::UndoMachine()
{
    events = std::list<Undoable *>();
    size = 0;
}

void UndoMachine::add_event(Undoable *event)
{
    events.push_back(event);
    size++;
}

void UndoMachine::insert_event(std::list<Undoable *>::iterator i, Undoable *event)
{
    std::list<Undoable *>::iterator j = i;
    j++;
    events.insert(j, event);
    size++;
}

std::list<Undoable *>::iterator UndoMachine::get_bookmark()
{
    std::list<Undoable *>::iterator bookmark = events.end();
    if (size > 0)
        bookmark--;
    else
        bookmark = events.begin();
    return bookmark;
}

void UndoMachine::undo()
{
    if (!events.empty())
    {
        Undoable *event = events.back();
        events.pop_back();
        event->undo();
        delete event;
        size--;
    }
}

void UndoMachine::undo(int num)
{
    while (num > 0)
        undo();
}

void UndoMachine::undo_all()
{
    undo_to(0);
}

void UndoMachine::undo_to(int to)
{
    while (size > to)
        undo();
}

void UndoMachine::clear_to(int to)
{
    while (size > to)
    {
        Undoable *event = events.back();
        events.pop_back();
        delete event;
        size--;
    }
    
}

int UndoMachine::num_events()
{
    return size;
}

void ContractEvent(UndoMachine *um, SPRNode *n, std::list<Undoable *>::iterator
                   bookmark)
{
    SPRNode *parent = n->parent();
    SPRNode *child;
    SPRNode *lc = n->lchild();
    SPRNode *rc = n->rchild();
    //SPRNode *ret = NULL;
    // contract out this SPRNode and give child to parent
    if (parent != NULL) {
        if (lc && !rc) {
            child = lc;
            um->add_event(new ChangeEdgePreInterval(child));
            um->add_event(new CutParent(child));
            um->add_event(new CutParent(n));
            if (n->is_protected() && !child->is_protected())
                um->add_event(new ProtectEdge(child));
        }
        else if (rc && !lc) {
            child = rc;
            um->add_event(new ChangeEdgePreInterval(child));
            um->add_event(new CutParent(child));
            um->add_event(new CutParent(n));
            if (n->is_protected() && !child->is_protected())
                um->add_event(new ProtectEdge(child));
        }
        else if (lc == NULL && rc == NULL) {
            um->insert_event(bookmark, new CutParent(n));
            parent->delete_child(n);
            ContractEvent(um, parent);
            parent->add_child(n);
        }
    }
    // if no parent then take children of single child and remove it
    else {
        
        // dead component or singleton, will be cleaned up by the forest
        if (n->get_children().empty()) {
            um->add_event(new ChangeName(n));
        }
        else if (n->get_children().size() == 1) {
            child = n->get_children().front();
            um->add_event(new CutParent(child));
            /* cluster hack - if we delete a cluster SPRNode then
             * we may try to use it later. This only happens once
             * per cluster so we can spend linear time to update
             * the forest
             */
            if (child->get_num_clustered_children() > 0) {
                //um->add_event(new CutParent(n));
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
                        um->add_event(new SetTwin(n));
                        um->add_event(new SetTwin(child->get_twin()));
                    }
                    um->add_event(new ChangeName(n));
                }
                um->add_event(new ChangePreNum(n));
                //um->add_event(new CutParent(n));
                std::list<SPRNode *>::iterator c;
                for(c = child->get_children().begin();
                    c != child->get_children().end();
                    c++) {
                    um->add_event(new CutParent(*c));
                }
                if (child->get_contracted_lc() != NULL)
                    um->add_event(new AddContractedLC(n));
                if (child->get_contracted_rc() != NULL)
                    um->add_event(new AddContractedRC(n));
            }
        }
    }
    
}

void ContractEvent(UndoMachine *um, SPRNode *n)
{
    std::list<Undoable *>::iterator bookmark = um->get_bookmark();
    ContractEvent(um, n, bookmark);
}




