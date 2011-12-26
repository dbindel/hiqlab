/*
 * matexpr-ir.h
 *   matexpr intermediate representation header.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef MATEXPR_IR_H
#define MATEXPR_IR_H

#include <string>
#include <list>

enum {
    IR_LABEL, IR_BINOP, IR_NEG, IR_VALUE, 
    IR_READ, IR_WRITE, IR_COPY, IR_CALL, IR_LINE_CPP
};

struct IRentry {
    IRentry(int tag) : 
        tag(tag), arg1(0), arg2(0), iscomplex(0) {
        static int slotcount = 0;
        slot = ++slotcount;
    }

    int tag;
    char op;
    std::string name;
    IRentry* arg1;
    IRentry* arg2;

    int slot;
    bool live;
    bool iscomplex;
};

typedef IRentry* IRop;
typedef std::list<IRop> IRlist;

void clear(IRlist& ir);
void fold_zeros(IRlist& ir);
void fold_identity(IRlist& ir);
void copy_forward(IRlist& ir);
void eliminate_common(IRlist& ir);
void eliminate_dead(IRlist& ir);
void relabel_ir(IRlist& ir);
void mark_complex(IRlist& ir);

#endif /* MATEXPR_IR_H */
