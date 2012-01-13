/*
 * matexpr-ast2const.cc
 *   Evaluate simple scalar constant expressions in AST.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "matexpr-ast.h"
#include <math.h>
#include <stdlib.h>

double ast2double(ASTNode* node, bool& is_const)
{
    if (!is_const)
        return 1;

    if (node->op == AST_ADD)
        return ast2double(node->l, is_const) + ast2double(node->r, is_const);
    else if (node->op == AST_SUB)
        return ast2double(node->l, is_const) - ast2double(node->r, is_const);
    else if (node->op == AST_MUL)
        return ast2double(node->l, is_const) * ast2double(node->r, is_const);
    else if (node->op == AST_DIV)
        return ast2double(node->l, is_const) / ast2double(node->r, is_const);
    else if (node->op == AST_NEG)
        return -ast2double(node->l, is_const);
    else if (node->op == AST_VALUE)
        return atof(node->name.c_str());
    else {
        is_const = false;
        return 1;
    }
}


int ast2int(ASTNode* node, bool& is_iconst)
{
    double x = ast2double(node, is_iconst);
    if (is_iconst) {
        int i = (int) x;
        if (x == (double) i)
            return i;
    }
    is_iconst = false;
    return 1;
}
