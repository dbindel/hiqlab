/*
 * matexpr-cwrite.cc
 *   matexpr code output implementation.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "matexpr-ir.h"

#include <set>
#include <algorithm>
#include <cassert>

using std::set;

extern "C" int get_linepos();


void clear(IRlist& ir)
{
    for (IRlist::iterator i = ir.begin(); i != ir.end(); ++i)
        delete(*i);
    ir.clear();
}


bool is_zero(IRop o)
{
    if (o->tag == IR_VALUE)
        return (atof(o->name.c_str()) == 0);
    if (o->tag == IR_COPY)
        return is_zero(o->arg1);
    if (o->tag == IR_NEG)
        return is_zero(o->arg1);
    if (o->tag == IR_BINOP)
        return
            (o->op == '+' && (is_zero(o->arg1) && is_zero(o->arg2))) ||
            (o->op == '-' && (is_zero(o->arg1) && is_zero(o->arg2))) ||
            (o->op == '*' && (is_zero(o->arg1) || is_zero(o->arg2))) ||
            (o->op == '/' && is_zero(o->arg1));
    return false;
}


bool is_one(IRop o)
{
    return
        (o->tag == IR_VALUE && atof(o->name.c_str()) == 1) ||
        (o->tag == IR_COPY && is_one(o->arg1));
}


void fold_zeros(IRlist& ir)
{
    for (IRlist::iterator i = ir.begin(); i != ir.end(); ++i) {
        IRop o = *i;
        if (is_zero(o)) {
            o->tag = IR_VALUE;
            o->name = "0";
        }
    }
}


void fold_identity(IRlist& ir)
{
    for (IRlist::iterator i = ir.begin(); i != ir.end(); ++i) {
        IRop o = *i;
        if (o->tag == IR_BINOP) {
            if (o->op == '*' && is_one(o->arg1) ||
                o->op == '+' && is_zero(o->arg1)) {
                o->tag = IR_COPY;
                o->arg1 = o->arg2;
                o->arg2 = NULL;
            } else if (o->op == '*' && is_one(o->arg2) ||
                       o->op == '+' && is_zero(o->arg2) ||
                       o->op == '-' && is_zero(o->arg2) ||
                       o->op == '/' && is_one(o->arg2)) {
                o->tag = IR_COPY;
            }
        }
    }
}


IRop copy_forward(IRop o)
{
    if (o->tag == IR_COPY)
        return copy_forward(o->arg1);
    return o;
}


void copy_forward(IRlist& ir)
{
    for (IRlist::iterator i = ir.begin(); i != ir.end(); ++i) {
        IRop o = *i;
        if (o->tag == IR_BINOP) {
            o->arg1 = copy_forward(o->arg1);
            o->arg2 = copy_forward(o->arg2);
        } else if (o->tag == IR_NEG || o->tag == IR_WRITE) {
            o->arg1 = copy_forward(o->arg1);
        }
    }
}


struct op_compare {
    bool operator()(const IRop op1, const IRop op2) const
    {
        // Same operation?
        if (op1->tag != op2->tag)
            return (op1->tag < op2->tag);
        if (op1->tag == IR_BINOP && op1->op != op2->op)
            return (op1->op < op2->op);
        
        // Same arguments?
        if (op1->tag == IR_VALUE) { 

            return (atof(op1->name.c_str()) < atof(op2->name.c_str()));

        } else if (op1->tag == IR_READ) {

            return (op1->name < op2->name);

        } else if (op1->tag == IR_BINOP && 
                   (op1->op == '+' || op1->op == '*')) {

            // Associative ops -- sort the slots before comparing
            int op1a = std::min(op1->arg1->slot, op1->arg2->slot);
            int op1b = std::max(op1->arg1->slot, op1->arg2->slot);
            int op2a = std::min(op2->arg1->slot, op2->arg2->slot);
            int op2b = std::max(op2->arg1->slot, op2->arg2->slot);
            if (op1a != op2a)
                return (op1a < op2a);
            if (op1b != op2b)
                return (op1b < op2b);
            return false;

        } else if (op1->tag == IR_BINOP) {

            if (op1->arg1->slot != op2->arg1->slot)
                return (op1->arg1->slot < op2->arg1->slot);
            if (op1->arg2->slot != op2->arg2->slot)
                return (op1->arg2->slot < op2->arg2->slot);
            return false;

        } else if (op1->tag == IR_NEG) {

            return (op1->arg1->slot < op2->arg1->slot);

        } else if (op1->tag == IR_CALL) {

            if (op1->name != op2->name)
                return (op1->name < op2->name);
            return (op1->arg1->slot < op2->arg1->slot);

        }
        return (op1->slot < op2->slot);
    }
};


void eliminate_common(IRlist& ir)
{
    set<IRop,op_compare> subexprs;
    for (IRlist::iterator i = ir.begin(); i != ir.end(); ++i) {

        IRop o = *i;
        if (o->tag == IR_BINOP) {
            o->arg1 = copy_forward(o->arg1);
            o->arg2 = copy_forward(o->arg2);
        } else if (o->tag == IR_NEG || o->tag == IR_WRITE) {
            o->arg1 = copy_forward(o->arg1);
        }

        set<IRop>::iterator subexpr = subexprs.find(o);
        if (subexpr != subexprs.end()) {
            o->tag = IR_COPY;
            o->arg1 = *subexpr;
        } else {
            subexprs.insert(o);
        }
    }
}


void mark_live(IRop o)
{
    o->live = true;
    if (o->tag == IR_BINOP) {
        mark_live(o->arg1);
        mark_live(o->arg2);
    } else if (o->tag == IR_NEG || o->tag == IR_COPY || 
               o->tag == IR_WRITE || o->tag == IR_CALL) {
        mark_live(o->arg1);
    }
}


void eliminate_dead(IRlist& ir)
{
   for (IRlist::iterator i = ir.begin(); i != ir.end(); ++i) {
        IRop o = *i;
        o->live = false;
        if (o->tag == IR_WRITE || o->tag == IR_LABEL || o->tag == IR_LINE_CPP)
            mark_live(o);
   }

   for (IRlist::iterator i = ir.begin(); i != ir.end();) {
       IRop o = *i;
       if (!o->live) {
           IRlist::iterator dead = i;
           ++i;
           ir.erase(dead);
           delete o;
       } else {
           ++i;
       }
   }
}


void relabel_ir(IRlist& ir)
{
    int slots = 0;
    for (IRlist::iterator i = ir.begin(); i != ir.end(); ++i)
        (*i)->slot = ++slots;
}


void mark_complex(IRlist& ir)
{
    for (IRlist::iterator i = ir.begin(); i != ir.end(); ++i) {
        IRop o = *i;
        if (o->tag == IR_CALL) {
            if (o->name == "real" || o->name == "imag" || o->name == "abs")
                o->iscomplex = false;
            else
                o->iscomplex = o->arg1->iscomplex;
        } else if (o->tag == IR_COPY || o->tag == IR_NEG)
            o->iscomplex = o->arg1->iscomplex;
        else if (o->tag == IR_BINOP)
            o->iscomplex = (o->arg1->iscomplex || o->arg2->iscomplex);
    }
}
