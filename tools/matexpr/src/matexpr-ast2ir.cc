/*
 * matexpr-ast2ir.cc
 *   matexpr code generation implementation.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "matexpr-ast2ir.h"
#include <cassert>
#include <cstdio>

#define ME IRWriter

using std::string;


void ME::operator()(ASTlist& stmts)
{
    // Fetch data and perform computations
    for (ASTlist::iterator i = stmts.begin(); i != stmts.end(); ++i) {
	if ((*i)->op == AST_ASSIGN) 
	    assign(*i);
	else if ((*i)->op == AST_INPUT ||
		 (*i)->op == AST_INOUT)
	    input(*i, false);
        else if ((*i)->op == AST_INPUTZ || 
                 (*i)->op == AST_INOUTZ)
            input(*i, true);
    }

    // Transfer data back to external environment
    for (ASTlist::iterator i = stmts.begin(); i != stmts.end(); ++i) {
	if ((*i)->op == AST_OUTPUT || 
	    (*i)->op == AST_INOUT ||
            (*i)->op == AST_INOUTZ) 
	    output(*i);
    }
}


IRop ME::add(IRop o)
{
    ir.push_back(o);
    return o;
}


void ME::assign(ASTNode* node)
{
    gen_label("/* Assignment at %d */\n", node->line);
    values[node->l->name] = expr(node->r);
}


void ME::input(ASTNode* node, bool iscomplex)
{
    gen_label("/* Input at %d */\n", node->line);
    if (node->array_flag == 1) {
	int m = node->m;
	int n = node->n;
	Matrix result(m,n);
	for (int j = 0; j < n; ++j)
	    for (int i = 0; i < m; ++i)
                result(i,j) = gen_read(node->l->name, node->lda, i, j,
                                       iscomplex);
	values[node->l->name] = result;
    } else if (node->array_flag == 2) {
	int m = node->m;
	Matrix result(m,m);
	for (int j = 0; j < m; ++j)
	    for (int i = 0; i <= j; ++i)
                result(i,j) = result(j,i) =
                    gen_read(node->l->name, node->lda, i, j, iscomplex);
	values[node->l->name] = result;
    } else {
	Matrix result(1,1);
	result(0,0) = gen_read(node->l->name, iscomplex);
	values[node->l->name] = result;
    }
}


void ME::output(ASTNode* node)
{
    gen_label("/* Output at %d */\n", node->line);
    Matrix result = values[node->l->name];
    if (node->array_flag) {
	for (int j = 0; j < node->n; ++j)
	    for (int i = 0; i < node->m; ++i)
                gen_write(node->l->name, node->lda, i, j, result(i,j));

    } else {
        gen_write(node->l->name, result(0,0));
    }
}


ME::Matrix ME::expr(ASTNode* node)
{
    if (!node)
	return Matrix();

    if      (node->op == AST_ID)         return id(node);
    else if (node->op == AST_VALUE)      return value(node);
    else if (node->op == AST_RANGE)      return range(node);
    else if (node->op == AST_ADD)        return binop(node, '+');
    else if (node->op == AST_SUB)        return binop(node, '-');
    else if (node->op == AST_MUL)        return mul(node);
    else if (node->op == AST_DIV)        return div(node);
    else if (node->op == AST_NEG)        return neg(node);
    else if (node->op == AST_TRANSP)     return transp(node);
    else if (node->op == AST_VCAT)       return vcat(node);
    else if (node->op == AST_HCAT)       return hcat(node);
    else if (node->op == AST_CALL)       return call(node);
    else if (node->op == AST_SUBSCRIPT)  return subscript(node);
    else
	assert(0);
    return Matrix();
}


ME::Matrix ME::id(ASTNode* node)
{
    return values[node->name];
}


ME::Matrix ME::value(ASTNode* node)
{
    Matrix result(1,1);
    result(0,0) = gen_value(node->name);
    return result;
}


ME::Matrix ME::binop(ASTNode* node, char opc)
{
    if (node->l->m == 1 && node->l->n == 1)
	return binop_scalar_mat(node, opc);
    else if (node->r->m == 1 && node->r->n == 1)
	return binop_mat_scalar(node, opc);

    Matrix xl = expr(node->l);
    Matrix xr = expr(node->r);
    Matrix result(xl.m, xl.n);
    for (int j = 0; j < result.n; ++j)
	for (int i = 0; i < result.m; ++i)
            result(i,j) = gen_binop(opc, xl(i,j), xr(i,j));
    return result;
}


ME::Matrix ME::binop_scalar_mat(ASTNode* node, char opc)
{
    Matrix xl = expr(node->l);
    Matrix xr = expr(node->r);
    Matrix result(xr.m, xr.n);
    for (int j = 0; j < result.n; ++j) 
	for (int i = 0; i < result.m; ++i) 
            result(i,j) = gen_binop(opc, xl(0,0), xr(i,j));
    return result;
}


ME::Matrix ME::binop_mat_scalar(ASTNode* node, char opc)
{
    Matrix xl = expr(node->l);
    Matrix xr = expr(node->r);
    Matrix result(xl.m, xl.n);
    for (int j = 0; j < result.n; ++j) 
	for (int i = 0; i < result.m; ++i) 
            result(i,j) = gen_binop(opc, xl(i,j), xr(0,0));
    return result;
}


ME::Matrix ME::mul(ASTNode* node)
{
    if (node->l->m == 1 && node->l->n == 1)
	return binop_scalar_mat(node, '*');
    else if (node->r->m == 1 && node->r->n == 1)
	return binop_mat_scalar(node, '*');

    Matrix xl = expr(node->l);
    Matrix xr = expr(node->r);
    Matrix result(xl.m, xr.n);
    for (int j = 0; j < result.n; ++j) {
	for (int i = 0; i < result.m; ++i) {
            result(i,j) = gen_binop('*', xl(i,0), xr(0,j));
	    for (int k = 1; k < xl.n; ++k)
                result(i,j) = gen_binop('+', result(i,j), 
                                        gen_binop('*', xl(i,k), xr(k,j)));
	}
    }
    return result;
}


ME::Matrix ME::div(ASTNode* node)
{
    return binop_mat_scalar(node, '/');
}


ME::Matrix ME::neg(ASTNode* node)
{
    Matrix xl = expr(node->l);
    Matrix result(xl.m, xl.n);
    for (int j = 0; j < result.n; ++j) 
	for (int i = 0; i < result.m; ++i) 
            result(i,j) = gen_neg(xl(i,j));
    return result;
}


ME::Matrix ME::range(ASTNode* node)
{
    bool is_iconst = true;
    int lo = ast2int(node->l, is_iconst);
    int hi = ast2int(node->r, is_iconst);
    int n = hi-lo+1;
    Matrix result(n,1);
    for (int ii = 0; ii < n; ++ii) {
        char buf[128];
        sprintf(buf, "%d", lo+ii);
        result(ii,0) = gen_value(buf);
    }
    return result;
}


ME::Matrix ME::transp(ASTNode* node)
{
    Matrix xl = expr(node->l);
    Matrix result(xl.n, xl.m);
    for (int j = 0; j < xl.n; ++j)
	for (int i = 0; i < xl.m; ++i)
	    result(j,i) = xl(i,j);
    return result;
}


ME::Matrix ME::vcat(ASTNode* node)
{
    Matrix xl = expr(node->l);
    Matrix xr = expr(node->r);
    Matrix result(xl.m + xr.m, xl.n);
    for (int j = 0; j < xl.n; ++j) {
	for (int i = 0; i < xl.m; ++i)
	    result(i,j) = xl(i,j);
	for (int i = 0; i < xr.m; ++i)
	    result(xl.m+i,j) = xr(i,j);
    }
    return result;
}


ME::Matrix ME::hcat(ASTNode* node)
{
    Matrix xl = expr(node->l);
    Matrix xr = expr(node->r);
    Matrix result(xl.m, xl.n + xr.n);
    for (int j = 0; j < xl.n; ++j)
	for (int i = 0; i < xl.m; ++i)
	    result(i,j) = xl(i,j);
    for (int j = 0; j < xr.n; ++j)
	for (int i = 0; i < xr.m; ++i)
	    result(i,xl.n+j) = xr(i,j);
    return result;
}


ME::Matrix ME::call(ASTNode* node)
{
    if (node->name == "eye")
        return eye(node);
    if (node->name == "deriv")
        return deriv(node);
    if (funtable[node->name])
        return macro_expand(node);

    Matrix arg = expr(node->l->l);
    Matrix result(arg.m, arg.n);
    for (int j = 0; j < arg.n; ++j)
	for (int i = 0; i < arg.m; ++i)
	    result(i,j) = gen_call(node->name, arg(i,j));
    return result;
}


ME::Matrix ME::subscript(ASTNode* node)
{
    Matrix result(1,1);
    Matrix array = values[node->name];
    ASTNode* args = node->l;
    bool is_iconst = true;
    int i = ast2int(args->l, is_iconst);
    int j = (args->r ? ast2int(args->r->l, is_iconst) : 1);
    result(0,0) = array(i-1,j-1);
    return result;
}


ME::Matrix ME::macro_expand(ASTNode* node)
{
    ME call_context(ir);

    // Bind formal arguments
    ASTNode* f = funtable[node->name];
    ASTNode* formals = f->l;
    ASTNode* actuals = node->l;
    while (formals && actuals) {
        call_context.values[formals->l->name] = expr(actuals->l);
        formals = formals->r;
        actuals = actuals->r;
    }

    // Emit code
    f->r->m = node->m;
    f->r->n = node->n;
    return call_context.expr(f->r);
}


ME::Matrix ME::eye(ASTNode* node)
{
    IRop one = gen_value("1");
    IRop zero = gen_value("0");
    Matrix result(node->m, node->m);
    for (int j = 0; j < node->m; ++j)
        for (int i = 0; i < node->m; ++i)
            result(i,j) = (i == j ? one : zero);
    return result;
}


ME::Matrix ME::deriv(ASTNode* node)
{
    Matrix result(node->m, node->n);
    Matrix fun  = expr(node->l->l);
    Matrix vars = expr(node->l->r->l);
    for (int jv = 0; jv < vars.n; ++jv)
        for (int iv = 0; iv < vars.m; ++iv)
            for (int jfn = 0; jfn < fun.n; ++jfn)
                for (int ifn = 0; ifn < fun.m; ++ifn)
                    result(iv*fun.m+ifn, jv*fun.n+jfn) = 
                        deriv(fun(ifn,jfn), vars(iv,jv));
    return result;
}


IRop ME::deriv(IRop f, IRop x)
{
    // f is something; x is a value
    if (f->tag == IR_BINOP && (f->op == '+' || f->op == '-')) {
        return gen_binop(f->op, deriv(f->arg1, x), deriv(f->arg2, x));
    } else if (f->tag == IR_NEG) {
        return gen_neg(deriv(f->arg1, x));
    } else if (f->tag == IR_BINOP && f->op == '*') {
        return gen_binop('+',
                         gen_binop('*', f->arg1, deriv(f->arg2, x)),
                         gen_binop('*', f->arg2, deriv(f->arg1, x)));
    } else if (f->tag == IR_BINOP && f->op == '/') {
        IRop num = 
            gen_binop('-',
                      gen_binop('*', f->arg2, deriv(f->arg1, x)),
                      gen_binop('*', f->arg1, deriv(f->arg2, x)));
        IRop denom =
            gen_binop('*', f->arg2, f->arg2);
        return gen_binop('/', num, denom);
    } else if (f->tag == IR_READ) {
        return gen_value((f->name == x->name) ? "1" : "0");
    } else if (f->tag == IR_VALUE) {
        return gen_value("0");
    } else if (f->tag == IR_CALL) {
        if (f->name == "sqrt") {
            return gen_binop('*', gen_binop('/', gen_value("0.5"), f),
                             deriv(f->arg1, x));
        } else if (f->name == "cos") {
            return gen_binop('*', gen_neg(gen_call("sin", f->arg1)),
                             deriv(f->arg1, x));
        } else if (f->name == "sin") {
            return gen_binop('*', gen_call("cos", f->arg1),
                             deriv(f->arg1, x));
        } else if (f->name == "log") {
            return gen_binop('/', deriv(f->arg1, x), f->arg1);
        } else if (f->name == "exp") {
            return gen_binop('*', f, deriv(f->arg1, x));
        } else {
            fprintf(stderr, "Don't know how to differentiate %s\n", 
                    f->name.c_str());
            return gen_value("0");
        }
    } else {
       fprintf(stderr, "Don't know how to differentiate: tag %d\n", f->tag);
       return gen_value("0");
    }
}


void ME::gen_label(const char* fmt, int line)
{
    if (gen_labels_flag) {
        IRop entry = new IRentry(IR_LABEL);
        sprintf(linebuf, fmt, line);
        entry->name = linebuf;
        add(entry);
    }
    if (gen_line_flag) {
        IRop entry = new IRentry(IR_LINE_CPP);
        sprintf(linebuf, "#line %d \"%s\"", line, infname.c_str());
        entry->name = linebuf;
        add(entry);
    }
}


IRop ME::gen_binop(char opc, IRop arg1, IRop arg2)
{
    IRop entry = new IRentry(IR_BINOP);
    entry->op = opc;
    entry->arg1 = arg1;
    entry->arg2 = arg2;
    return add(entry);
}


IRop ME::gen_neg(IRop arg)
{
    IRop entry = new IRentry(IR_NEG);
    entry->op = '-';
    entry->arg1 = arg;
    return add(entry);
}


IRop ME::gen_call(const string& name, IRop arg)
{
    IRop entry = new IRentry(IR_CALL);
    entry->name = name;
    entry->arg1 = arg;
    return add(entry);
}


IRop ME::gen_value(const string& name)
{
    IRop entry = new IRentry(IR_VALUE);
    entry->name = name;
    return add(entry);
}


IRop ME::gen_read(const string& name, bool iscomplex)
{
    IRop entry = new IRentry(IR_READ);
    entry->name = name;
    entry->iscomplex = iscomplex;
    return add(entry);
}


IRop ME::gen_read(const string& name, const string& lda, int i, int j, 
                  bool iscomplex)
{
    IRop entry = new IRentry(IR_READ);
    entry->name = array_location(name, lda, i, j);
    entry->iscomplex = iscomplex;
    return add(entry);
}


void ME::gen_write(const string& name, IRop result)
{
    IRop entry = new IRentry(IR_WRITE);
    entry->name = name;
    entry->arg1 = result;
    add(entry);
}


void ME::gen_write(const string& name, const string& lda, int i, int j, 
                   IRop result)
{
    IRop entry = new IRentry(IR_WRITE);
    entry->name = array_location(name, lda, i, j);
    entry->arg1 = result;
    add(entry);
}


string ME::array_location(const string& name, const string& lda, int i, int j)
{
    sprintf(linebuf, "%s[%d*%s+%d]", name.c_str(), j, lda.c_str(), i);
    return linebuf;
}
