/*
 * matexpr-ir2c.h
 *   Convert matexpr IR to C code.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef MATEXPR_IR2C_H
#define MATEXPR_IR2C_H

#include "matexpr-ir.h"
#include <stdio.h>


class CWriter {
public:
    CWriter(FILE* fp) : fp(fp) {}
    void operator()(IRlist& ir);

private:
    FILE* fp;
    void label   (IRop o);
    void binop   (IRop o);
    void neg     (IRop o);
    void value   (IRop o);
    void read    (IRop o);
    void write   (IRop o);
    void copy    (IRop o);
    void call    (IRop o);
    void line_cpp(IRop o);

    void output_assign(IRop o);
    void output_result(IRop o);
    void tab_line();
};

#endif /* MATEXPR_IR2C_H */
