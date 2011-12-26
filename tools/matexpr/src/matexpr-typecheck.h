/*
 * matexpr-typecheck.h
 *   matexpr type checker header.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef MATEXPR_TYPECHECKER_H
#define MATEXPR_TYPECHECKER_H

#include "matexpr-ast.h"

#include <stdio.h>
#include <string>
#include <map>


class TypeChecker {
public:
    TypeChecker(FILE* outfp) : outfp(outfp), errs(0) {}
    int operator()(ASTlist& stmts);

    FILE* get_outfp()  { return outfp; }

    void report_error() {
        ++errs;
        fprintf(outfp, "%d: ", line);
    }

private:
    FILE* outfp;
    int errs;
    int line;
    std::map<std::string, ASTNode*> values;

    void assign(ASTNode* node);
    void input(ASTNode* node);
    void output(ASTNode* node);
    
    ASTNode* expr(ASTNode* node);
    void id(ASTNode* node);
    void value(ASTNode* node);
    void binop(ASTNode* node, char op);
    void mul(ASTNode* node);
    void div(ASTNode* node);
    void neg(ASTNode* node);
    void range(ASTNode* node);
    void transp(ASTNode* node);
    void vcat(ASTNode* node);
    void hcat(ASTNode* node);
    void call(ASTNode* node);
    void subscript(ASTNode* node);
    void macro_expand(ASTNode* node);
    void eye(ASTNode* node);
    void deriv(ASTNode* node);
};

#endif /* MATEXPR_TYPECHECKER_H */
