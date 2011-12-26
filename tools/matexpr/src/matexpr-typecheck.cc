/*
 * matexpr-typecheck.cc
 *   matexpr type checker implementation.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "matexpr-typecheck.h"

#define ME TypeChecker


int ME::operator()(ASTlist& stmts)
{
    // Fetch data and perform computations
    for (ASTlist::iterator i = stmts.begin(); i != stmts.end(); ++i) {
	if ((*i)->op == AST_ASSIGN) 
	    assign(*i);
	else if ((*i)->op == AST_INPUT || (*i)->op == AST_INPUTZ ||
		 (*i)->op == AST_INOUT || (*i)->op == AST_INOUTZ)
	    input(*i);
    }

    // Transfer data back to external environment
    for (ASTlist::iterator i = stmts.begin(); i != stmts.end(); ++i) {
	if ((*i)->op == AST_OUTPUT || 
	    (*i)->op == AST_INOUT || (*i)->op == AST_INOUTZ) 
	    output(*i);
    }

    return errs;
}


void ME::assign(ASTNode* node)
{
    line = node->line;
    values[node->l->name] = expr(node->r);
    node->m = node->r->m;
    node->n = node->r->n;
}


void ME::input(ASTNode* node)
{
    line = node->line;
    values[node->l->name] = node;
}


void ME::output(ASTNode* node)
{
    line = node->line;
    ASTNode* result = values[node->l->name];
    if (!result) {
        report_error();
	fprintf(outfp, "%s not bound\n", node->l->name.c_str());
    } else if (node->m != result->m || node->n != result->n) {
        report_error();
	fprintf(outfp, "Size error (expected %d-by-%d, saw %d-by-%d)\n", 
		node->m, node->n, result->m, result->n);
    }
}


ASTNode* ME::expr(ASTNode* node)
{
    if (!node)
	return NULL;

    if      (node->op == AST_ID)      id(node);
    else if (node->op == AST_VALUE)   value(node);
    else if (node->op == AST_ADD)     binop(node, '+');
    else if (node->op == AST_SUB)     binop(node, '-');
    else if (node->op == AST_MUL)     mul(node);
    else if (node->op == AST_DIV)     div(node);
    else if (node->op == AST_NEG)     neg(node);
    else if (node->op == AST_RANGE)   range(node);
    else if (node->op == AST_TRANSP)  transp(node);
    else if (node->op == AST_VCAT)    vcat(node);
    else if (node->op == AST_HCAT)    hcat(node);
    else if (node->op == AST_CALL)    call(node);
    else if (node->op == AST_SUBSCRIPT) subscript(node);
    else {
	fprintf(outfp, "%d (internal): node type %d not yet implemented\n", 
		line, node->op);
	++errs;
	node->m = 0;
	node->n = 0;
    }
    return node;
}


void ME::id(ASTNode* node)
{
    ASTNode* value = values[node->name];
    if (value) {
	node->m = value->m;
	node->n = value->n;
    } else {
        report_error();
	fprintf(outfp, "%s not bound\n", node->name.c_str());
	node->m = 0;
	node->n = 0;
    }
}


void ME::value(ASTNode* node)
{
    node->m = 1;
    node->n = 1;
}


void ME::binop(ASTNode* node, char opc)
{
    ASTNode* xl = expr(node->l);
    ASTNode* xr = expr(node->r);
    if (xl->m == 1 && xl->n == 1) {
	node->m = xr->m;
	node->n = xr->n;
    } else if (xr->m == 1 && xr->n == 1) {
	node->m = xl->m;
	node->n = xl->n;
    } else if (xl->m != xr->m || xl->n != xr->n) {
        report_error();
	fprintf(outfp, "Type error: %d-by-%d %c %d-by-%d\n",
		xl->m, xl->n, opc, xr->m, xr->n);
	node->m = xl->m;
	node->n = xl->n;
    } else {
	node->m = xl->m;
	node->n = xl->n;
    }
}


void ME::mul(ASTNode* node)
{
    ASTNode* xl = expr(node->l);
    ASTNode* xr = expr(node->r);
    if (xl->m == 1 && xl->n == 1) {
	node->m = xr->m;
	node->n = xr->n;
    } else if (xr->m == 1 && xr->n == 1) {
	node->m = xl->m;
	node->n = xl->n;
    } else if (xl->n != xr->m) {
        report_error();
	fprintf(outfp, "Type error: %d-by-%d * %d-by-%d\n",
		xl->m, xl->n, xr->m, xr->n);
	node->m = xl->m;
	node->n = xr->n;
    } else {
	node->m = xl->m;
	node->n = xr->n;
    }
}


void ME::div(ASTNode* node)
{
    ASTNode* xl = expr(node->l);
    ASTNode* xr = expr(node->r);
    node->m = xl->m;
    node->n = xl->n;
    if (xr->m != 1 || xr->n != 1) {
        report_error();
	fprintf(outfp, "Type error: %d-by-%d / %d-by-%d\n",
		xl->m, xl->n, xr->m, xr->n);
    }
}


void ME::neg(ASTNode* node)
{
    ASTNode* xl = expr(node->l);
    node->m = xl->m;
    node->n = xl->n;
}


void ME::range(ASTNode* node)
{
    bool is_iconst = true;
    int lo = ast2int(node->l, is_iconst);
    int hi = ast2int(node->r, is_iconst);
    node->m = 0;
    node->n = 1;
    if (!is_iconst) {
        report_error();
        fprintf(outfp, "Range bounds must be integer constants\n");
    } else if (hi < lo) {
        report_error();
        fprintf(outfp, "Lower bound cannot be above upper bound\n");
    } else {
        node->m = hi-lo+1;
    }
}


void ME::transp(ASTNode* node)
{
    ASTNode* xl = expr(node->l);
    node->m = xl->n;
    node->n = xl->m;
}


void ME::vcat(ASTNode* node)
{
    ASTNode* xl = expr(node->l);
    ASTNode* xr = expr(node->r);
    if (xl->n != xr->n) {
        report_error();
	fprintf(outfp, "Type error vcat(%d-by-%d, %d-by-%d)\n", 
                xl->m, xl->n, xr->m, xr->n);
    }
    node->m = xl->m + xr->m;
    node->n = xl->n;
}


void ME::hcat(ASTNode* node)
{
    ASTNode* xl = expr(node->l);
    ASTNode* xr = expr(node->r);
    if (xl->m != xr->m) {
        report_error();
	fprintf(outfp, "Type error hcat(%d-by-%d, %d-by-%d)\n", 
                xl->m, xl->n, xr->m, xr->n);
    }
    node->m = xl->m;
    node->n = xl->n + xr->n;
}


void ME::call(ASTNode* node)
{
    if (values[node->name]) {
        node->op = AST_SUBSCRIPT;
        subscript(node);
        return;
    }

    if (funtable[node->name]) {
        macro_expand(node);
        return;
    }

    if (node->name == "eye") {
        eye(node);
        return;
    }

    if (node->name == "deriv") {
        deriv(node);
        return;
    }

    if (!node->l || node->l->r) {
        report_error();
        fprintf(outfp, "Currently only handle calls with exactly one arg\n");
        node->m = 0;
        node->n = 0;
    } else {
        ASTNode* arg = expr(node->l->l);
        node->m = arg->m;
        node->n = arg->n;
    }
}


void ME::subscript(ASTNode* node)
{
    ASTNode* array = values[node->name];
    ASTNode* args = node->l;

    node->m = 1;
    node->n = 1;

    if (!args) {
        report_error();
        fprintf(outfp, "Cannot subscript %s with no indices\n", 
                node->name.c_str());
    } else if (!args->r) {
        ASTNode* arg1 = expr(args->l);

        bool is_iconst = true;
        int i = ast2int(arg1, is_iconst);

        if (!is_iconst) {
            report_error();
            fprintf(outfp, "Indices must be scalar constants\n");
        } else if (i < 1 || i > array->m * array->n) {
            report_error();
            fprintf(outfp, "Index out of range\n");
        } 
    } else {
        ASTNode* arg1 = expr(args->l);
        ASTNode* arg2 = expr(args->r->l);

        bool is_iconst = true;
        int i = ast2int(arg1, is_iconst);
        int j = ast2int(arg2, is_iconst);

        if (args->r->r) {
            report_error();
            fprintf(outfp, "Too many indices!\n");
        } else if (!is_iconst) {
            report_error();
            fprintf(outfp, "Indices must be scalar constants\n");
        } else if (i < 1 || i > array->m || j < 1 || j > array->n) {
            report_error();
            fprintf(outfp, "Index out of range\n");
        } 
    }
}


void ME::macro_expand(ASTNode* node)
{
    node->m = 0;
    node->n = 0;
    ME call_context(outfp);
    call_context.line = line;

    // Check and bind formal arguments
    ASTNode* f = funtable[node->name];
    ASTNode* formals = f->l;
    ASTNode* actuals = node->l;
    while (formals && actuals) {
        call_context.values[formals->l->name] = expr(actuals->l);
        formals = formals->r;
        actuals = actuals->r;
    }

    // Check counts
    if (formals) {
        report_error();
        fprintf(outfp, "Not enough arguments in call to %s\n", 
                node->name.c_str());
    }
    if (actuals) {
        report_error();
        fprintf(outfp, "Too many arguments in call to %s\n", 
                node->name.c_str());
    }

    // Do typecheck
    f->r->m = 0;
    f->r->n = 0;
    call_context.expr(f->r);
    node->m = f->r->m;
    node->n = f->r->n;
    errs += call_context.errs;
}


void ME::eye(ASTNode* node)
{
    ASTNode* args = node->l;

    node->m = 1;
    node->n = 1;

    if (!args) {
        report_error();
        fprintf(outfp, "eye must take an argument\n");
    } else if (!args->r) {
        bool is_iconst = true;
        int m = ast2int(args->l, is_iconst);
        if (!is_iconst) {
            report_error();
            fprintf(outfp, "Size must be scalar constants\n");
        } else if (m <= 0) {
            report_error();
            fprintf(outfp, "Size out of range\n");
        } 
        node->m = m;
        node->n = m;
    } else {
        report_error();
        fprintf(outfp, "Too many arguments to eye!\n");
    }
}


void ME::deriv(ASTNode* node)
{
    ASTNode* args = node->l;

    node->m = 1;
    node->n = 1;

    if (!args || !args->r) {
        report_error();
        fprintf(outfp, "deriv must take two arguments\n");
    } else if (!args->r->r) {
        // FIXME: Would be good to check that this is good for diff, too.
        expr(args->l);
        expr(args->r->l);
        node->m = args->l->m * args->r->l->m;
        node->n = args->l->n * args->r->l->n;
    } else {
        report_error();
        fprintf(outfp, "Too many arguments to eye!\n");
    }    
}
