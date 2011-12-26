/*
 * matexpr-ast.h
 *   Abstract syntax for matrix expressions.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef MATEXPR_AST_H
#define MATEXPR_AST_H

#include <string>
#include <list>
#include <map>


enum {
    AST_RANGE, AST_ADD, AST_SUB, AST_MUL, AST_DIV, AST_NEG, 
    AST_TRANSP, AST_VCAT, AST_HCAT, AST_ASSIGN, AST_FUNCTION,
    AST_INPUT, AST_OUTPUT, AST_INOUT, AST_INPUTZ, AST_INOUTZ, 
    AST_DIM, AST_CALL, AST_SUBSCRIPT, AST_ARG, AST_ID, AST_VALUE
};

struct ASTNode {

    ASTNode(const char* id) : 
        op(AST_ID), name(id), l(0), r(0) { }
    ASTNode(int op, ASTNode* l, ASTNode* r) : 
	op(op), l(l), r(r) { }
    ASTNode(int op, ASTNode* l) : 
	op(op), l(l), r(0) { }
    ASTNode(int op, const char* id) :
        op(op), name(id), l(0), r(0) { }
    ASTNode(const char* id, ASTNode* arg) :
        op(AST_CALL), name(id), l(arg), r(0) {}

    ~ASTNode() {
	if (l) delete l;
	if (r) delete r;
    }

    int op;
    std::string name;
    ASTNode* l;
    ASTNode* r;
    int line;

    int m, n;
    int array_flag; // 1 for array decl, 2 for symm array decl
    std::string lda;

private:
    ASTNode(const ASTNode&);
    ASTNode& operator=(const ASTNode&);
};

typedef std::list<ASTNode*> ASTlist;
typedef std::map<std::string,ASTNode*> ASTdict;

extern std::string infname;
extern int nogen_flag;
extern int gen_labels_flag;
extern int gen_line_flag;
extern int err_count;
extern ASTlist stmts;
extern ASTdict funtable;

double ast2double(ASTNode* node, bool& is_const);
int    ast2int(ASTNode* node, bool& is_iconst);

#endif /* MATEXPR_AST_H */
