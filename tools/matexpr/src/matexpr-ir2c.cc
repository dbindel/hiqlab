/*
 * matexpr-ir2c.cc
 *   Convert matexpr IR to C code.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "matexpr-ir2c.h"

#include <cassert>
extern "C" int get_linepos();
extern std::string complexname;

#define ME CWriter


void ME::operator()(IRlist& ir)
{
    tab_line();
    fprintf(fp, "/* <generated matexpr> */ {\n");

    for (IRlist::iterator i = ir.begin(); i != ir.end(); ++i) {
        IRop op = *i;
        switch (op->tag) {
        case IR_LABEL:     label(op);    break;
        case IR_BINOP:     binop(op);    break;
        case IR_NEG:       neg(op);      break;
        case IR_VALUE:     value(op);    break;
        case IR_READ:      read(op);     break;
        case IR_WRITE:     write(op);    break;
        case IR_COPY:      copy(op);     break;
        case IR_CALL:      call(op);     break;
        case IR_LINE_CPP:  line_cpp(op); break;
        default:        assert(0); break;
        }
    }

    tab_line();
    fprintf(fp, "} /* </generated matexpr> */\n");    
}


void ME::label(IRop o)
{
    fprintf(fp, "\n");
    tab_line();
    fprintf(fp, "%s", o->name.c_str());
}


void ME::binop(IRop o)
{
    output_assign(o);
    output_result(o->arg1);
    fprintf(fp, " %c ", o->op);
    output_result(o->arg2);
    fprintf(fp, ";\n");
}


void ME::neg(IRop o)
{
    output_assign(o);
    fprintf(fp, "-");
    output_result(o->arg1);
    fprintf(fp, ";\n");
}


void ME::value(IRop o)
{
/* Doesn't actually need to write anything because of constant folding */
/*
    output_assign(o);
    fprintf(fp, "%s;\n", o->name.c_str());
*/
}


void ME::read(IRop o)
{
    output_assign(o);
    fprintf(fp, "%s;\n", o->name.c_str());
}


void ME::write(IRop o)
{
    tab_line();
    fprintf(fp, o->name.c_str());
    fprintf(fp, " = ");
    output_result(o->arg1);
    fprintf(fp, ";\n");
}


void ME::copy(IRop o)
{
    output_assign(o);
    output_result(o->arg1);
    fprintf(fp, ";\n");
}


void ME::call(IRop o)
{
    output_assign(o);
    fprintf(fp, "%s(", o->name.c_str());
    output_result(o->arg1);
    fprintf(fp, ");\n");
}


void ME::line_cpp(IRop o)
{
    fprintf(fp, "%s\n", o->name.c_str());
}


void ME::output_assign(IRop o)
{
    tab_line();
    if (o->iscomplex)
        fprintf(fp, "%s tmp%d_ = ", complexname.c_str(), o->slot);
    else
        fprintf(fp, "double tmp%d_ = ", o->slot);
}


void ME::output_result(IRop o)
{
    if (o->tag == IR_VALUE)
        fprintf(fp, "%s", o->name.c_str());
    else
        fprintf(fp, "tmp%d_", o->slot);
}


void ME::tab_line()
{
    for (int j = 0; j < get_linepos(); ++j)
        fputc(' ', fp);
}
