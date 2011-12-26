/*
 * matexpr-ast2ir.h
 *   matexpr code generation header.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#ifndef MATEXPR_CGEN_H
#define MATEXPR_CGEN_H

#include "matexpr-ast.h"
#include "matexpr-ir.h"

#include <string>
#include <vector>
#include <map>


class IRWriter {
public:
    IRWriter(IRlist& ir) : ir(ir) {}
    void operator()(ASTlist& stmts);

    struct Matrix {
	Matrix() : m(0), n(0), data(1) {}
	Matrix(int m, int n) : m(m), n(n), data(m*n) {}
	
	IRop& operator()(int i, int j) { return data[j*m+i]; }
	int m, n;
	std::vector<IRop> data;
    };

    IRop add(IRop o);

private:
    char linebuf[1024];
    std::map<std::string, Matrix> values;
    IRlist& ir;

    void assign(ASTNode* node);
    void input(ASTNode* node, bool iscomplex);
    void output(ASTNode* node);

    Matrix expr(ASTNode* node);
    Matrix id(ASTNode* node);
    Matrix value(ASTNode* node);
    Matrix binop(ASTNode* node, char op);
    Matrix binop_mat_scalar(ASTNode* node, char op);
    Matrix binop_scalar_mat(ASTNode* node, char op);
    Matrix mul(ASTNode* node);
    Matrix div(ASTNode* node);
    Matrix neg(ASTNode* node);
    Matrix range(ASTNode* node);
    Matrix transp(ASTNode* node);
    Matrix vcat(ASTNode* node);
    Matrix hcat(ASTNode* node);
    Matrix call(ASTNode* node);
    Matrix macro_expand(ASTNode* node);
    Matrix subscript(ASTNode* node);
    Matrix eye(ASTNode* node);
    Matrix deriv(ASTNode* node);
    IRop deriv(IRop f, IRop x);

    void gen_label(const char* fmt, int line);
    IRop gen_binop(char opc, IRop arg1, IRop arg2);
    IRop gen_neg  (IRop arg);
    IRop gen_call (const std::string& name, IRop arg);
    IRop gen_value(const std::string& name);
    IRop gen_read (const std::string& name, bool iscomplex);
    IRop gen_read (const std::string& name, const std::string& lda, 
                   int i, int j, bool iscomplex);
    void gen_write(const std::string& name, IRop result);
    void gen_write(const std::string& name, const std::string& lda, 
                   int i, int j, IRop result);

    std::string array_location(const std::string& name, const std::string& lda,
                               int i, int j);
};


#endif /* MATEXPR_CGEN_H */
