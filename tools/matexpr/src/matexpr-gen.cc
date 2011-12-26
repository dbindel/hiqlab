/*
 * matexpr-ast.cc
 *   Utilities for matexpr.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "matexpr-ast.h"
#include "matexpr-typecheck.h"
#include "matexpr-ast2ir.h"
#include "matexpr-ir.h"
#include "matexpr-ir2c.h"


void clear_stmts()
{
    for (ASTlist::iterator i = stmts.begin(); i != stmts.end(); ++i)
        delete(*i);
    stmts.clear();
}


extern "C"
void generate(FILE* outfp)
{
    TypeChecker checker(stderr);
    int local_errs = checker(stmts);
    err_count += local_errs;
    if (!nogen_flag && local_errs == 0 && outfp != NULL) {

        IRlist ir;
	IRWriter codegen(ir);
        codegen(stmts);

        fold_zeros(ir);
        fold_identity(ir);
        copy_forward(ir);
        eliminate_common(ir);
        eliminate_dead(ir);
        relabel_ir(ir);
        mark_complex(ir);

        CWriter cwriter(outfp);
        cwriter(ir);
        clear(ir);

    }
    clear_stmts();
}

