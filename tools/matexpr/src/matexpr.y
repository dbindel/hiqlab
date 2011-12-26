%{
/*
 * matexpr.y
 *   Parser for MATLAB-like expressions.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "matexpr-ast.h"

#include <stdlib.h>
#include <string.h>

#include <string>
#include <map>

using std::string;
using std::map;

extern "C" {
    int yylex();
    int yywrap();
    int yyerror(char* s);
}

int linenum;
int err_count;
int nogen_flag;
int gen_labels_flag;
int gen_line_flag;
int typecheck_flag;
static int which_iospec;

string infname;
string complexname;
string default_complexname;
ASTlist stmts;
ASTdict funtable;

%}

%union {
    struct ASTNode* ast;
    char* s;
}

%token INPUT OUTPUT INOUT COMPLEX SYMMETRIC
%token ADDTO PRIME GENBLOCK FUNCTION
%token <s> ID STRING
%token <s> INUMBER FNUMBER
%type <ast> expr expr2 matrix matdata rowdata 
%type <ast> rowrest colrest args argsrest formals formalsrest
%type <s> lda declbase

%left ':'
%left '-' '+'
%left '*' '/'
%nonassoc UMINUS 

%%
statements: statement statements | ;

statement: 
  GENBLOCK {
      complexname = default_complexname;
    } genparams '>' ;
  | ID '=' expr ';'   { 
      ASTNode* stmt = new ASTNode(AST_ASSIGN, new ASTNode($1), $3); 
      stmt->line = linenum;
      stmts.push_back(stmt);
      free($1);
    }
  | ID ADDTO expr ';' {
      ASTNode* stmt = new ASTNode(AST_ASSIGN, new ASTNode($1),
				  new ASTNode(AST_ADD, new ASTNode($1), $3));
      stmt->line = linenum;
      stmts.push_back(stmt);
      free($1);
    }
  | FUNCTION ID '(' formals ')' '=' expr ';' {
      ASTNode* fn = new ASTNode(AST_FUNCTION, $4, $7);
      funtable[$2] = fn; // FIXME: Clean up old defs
      free($2);
    }
  | iospec decllist ';' ;

genparams: 
    COMPLEX '=' STRING { 
      complexname = $3; 
      free($3); 
    }
  | ;

iospec: 
    INPUT          { which_iospec = AST_INPUT;  }
  | OUTPUT         { which_iospec = AST_OUTPUT; }
  | INOUT          { which_iospec = AST_INOUT;  }
  | COMPLEX INPUT  { which_iospec = AST_INPUTZ;  }
  | COMPLEX INOUT  { which_iospec = AST_INOUTZ;  }

decllist:
    decllist ',' decl 
  | decl ;

decl: 
    declbase { free($1); }
  | declbase '=' expr {
      ASTNode* stmt = new ASTNode(AST_ASSIGN, new ASTNode($1), $3); 
      stmt->line = linenum;
      stmts.push_back(stmt);
      free($1);
    }
  | declbase ADDTO expr {
      ASTNode* stmt = new ASTNode(AST_ASSIGN, new ASTNode($1), 
                                  new ASTNode(AST_ADD, new ASTNode($1), $3)); 
      stmt->line = linenum;
      stmts.push_back(stmt);
      free($1);
    } ;

declbase:
    ID { 
      char ldabuf[16];
      ASTNode* stmt = new ASTNode(which_iospec, new ASTNode($1)); 
      stmt->line = linenum;
      stmt->m = 1;
      stmt->n = 1;
      stmt->lda = "1";
      stmt->array_flag = 0;
      stmts.push_back(stmt);
  }
  | ID '(' INUMBER ')' { 
      ASTNode* stmt = new ASTNode(which_iospec, new ASTNode($1));
      stmt->line = linenum;
      stmt->m = atoi($3);
      stmt->n = 1;
      stmt->lda = $3;
      stmt->array_flag = 1;
      stmts.push_back(stmt);
      free($3);
  }
  | ID '(' INUMBER ',' INUMBER ')' {
      ASTNode* stmt = new ASTNode(which_iospec, new ASTNode($1));
      stmt->line = linenum;
      stmt->m = atoi($3);
      stmt->n = atoi($5);
      stmt->lda = $3;
      stmt->array_flag = 1;
      stmts.push_back(stmt);
      free($3);
      free($5);
  } 
  | ID SYMMETRIC '(' INUMBER ')' {
      ASTNode* stmt = new ASTNode(which_iospec, new ASTNode($1));
      stmt->line = linenum;
      stmt->m = stmt->n = atoi($4);
      stmt->lda = $4;
      stmt->array_flag = 2;
      stmts.push_back(stmt);
      free($4);
  }
  | ID '[' lda ']' '(' INUMBER ',' INUMBER ')' {
      ASTNode* stmt = new ASTNode(which_iospec, new ASTNode($1));
      stmt->line = linenum;
      stmt->m = atoi($6);
      stmt->n = atoi($8);
      stmt->lda = $3;
      stmt->array_flag = 1;
      stmts.push_back(stmt);
      free($3);
      free($6);
      free($8);
  } ;

lda:
    ID      { $$ = $1; }
  | INUMBER { $$ = $1; } ;

expr:
    expr ':' expr { $$ = new ASTNode(AST_RANGE, $1, $3); }
  | expr '+' expr { $$ = new ASTNode(AST_ADD, $1, $3); }
  | expr '-' expr { $$ = new ASTNode(AST_SUB, $1, $3); }
  | expr '*' expr { $$ = new ASTNode(AST_MUL, $1, $3); }
  | expr '/' expr { $$ = new ASTNode(AST_DIV, $1, $3); }
  | '-' expr  %prec UMINUS 
                  { $$ = new ASTNode(AST_NEG, $2); } 
  | expr2 ;

expr2:
    expr2 PRIME   { $$ = new ASTNode(AST_TRANSP, $1); }
  | '(' expr ')'  { $$ = $2; }
  | ID            { $$ = new ASTNode(AST_ID,    $1); free($1); }
  | INUMBER       { $$ = new ASTNode(AST_VALUE, $1); free($1); 
                    $$->name += ".0"; } 
  | FNUMBER       { $$ = new ASTNode(AST_VALUE, $1); free($1); } 
  | matrix        { $$ = $1; } 
  | ID '(' args ')' { 
      $$ = new ASTNode($1, $3);
      free($1);
  } ;

matrix: 
    '[' matdata ']' { $$ = $2; } ;

matdata: 
    rowdata rowrest { $$ = new ASTNode(AST_VCAT, $1, $2); }
  | rowdata         { $$ = $1;   }
  |                 { $$ = NULL; } ;

rowrest:
    ';' rowdata rowrest { $$ = new ASTNode(AST_VCAT, $2, $3); }
  | ';' rowdata         { $$ = $2; } ;

rowdata:
    expr colrest     { $$ = new ASTNode(AST_HCAT, $1, $2); } ;
  | expr ;

colrest:
    ',' expr colrest { $$ = new ASTNode(AST_HCAT, $2, $3); }
  | ',' expr         { $$ = $2; } ;

args:
    expr argsrest { $$ = new ASTNode(AST_ARG, $1, $2); }
  | expr          { $$ = new ASTNode(AST_ARG, $1);     }
  |               { $$ = NULL; } ;

argsrest:
    ',' expr argsrest { $$ = new ASTNode(AST_ARG, $2, $3);   }
  | ',' expr          { $$ = new ASTNode(AST_ARG, $2);       } ;

formals:
    ID formalsrest { $$ = new ASTNode(AST_ARG, new ASTNode($1), $2); }
  | ID             { $$ = new ASTNode(AST_ARG, new ASTNode($1));     }
  |                { $$ = NULL; } ;

formalsrest:
    ',' ID formalsrest { $$ = new ASTNode(AST_ARG, new ASTNode($2), $3); }
  | ',' ID             { $$ = new ASTNode(AST_ARG, new ASTNode($2));     } ;

%%
#include <stdio.h>
#include <string.h>

extern FILE* yyin;
FILE* outfp;
extern "C" void start_filepos();
extern "C" void print_cpp_line(int linenum);

int yywrap()
{
    return 1;
}

int yyerror(char* s)
{
    fprintf(stderr, "Parse error (%s: %d): %s\n", 
	    infname.c_str(), linenum, s);
}

void print_cpp_line(int linenum)
{
    if (outfp && gen_line_flag)
        fprintf(outfp, "#line %d \"%s\"\n", linenum, infname.c_str());
}

int main(int argc, char** argv)
{
    outfp = stdout;
    nogen_flag = 0;
    err_count = 0;
    gen_labels_flag = 0;
    gen_line_flag = 0;
    default_complexname = "std::complex<double>";
    if (argc == 1) {
	infname = "<stdin>";
        start_filepos();
	yyparse();
    } else {
	for (int i = 1; i < argc; ++i) {
	    if (strcmp(argv[i], "-nogen") == 0) {
		nogen_flag = 1;
	    } else if (strcmp(argv[i], "-check") == 0) {
		typecheck_flag = 1;
		outfp = NULL;
            } else if (strcmp(argv[i], "-comment") == 0) {
                gen_labels_flag = 1;
            } else if (strcmp(argv[i], "-line") == 0) {
                gen_line_flag = 1;
            } else if (strcmp(argv[i], "-c99complex") == 0) {
                default_complexname = "double _Complex";
	    } else {
		yyin = fopen(argv[i], "r");
		if (!yyin) {
		    fprintf(stderr, "Could not open %s for input\n", argv[i]);
		    exit(-1);
 		}
		infname = argv[i];
                start_filepos();
		yyparse();
		fclose(yyin);
	    }
	}
    }
    return err_count;
}
