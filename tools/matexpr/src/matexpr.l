%{
/*
 * matexpr.l
 *   Lexer for MATLAB-like expressions.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */

#include "matexpr.hh"
#include <string.h>

void generate(FILE* outfp);
void print_cpp_line(int linenum);
extern int linenum;
static int linepos;
extern FILE* outfp;

int  get_linepos()   { return linepos; }
void start_filepos() { linenum = 1; linepos = 0; print_cpp_line(linenum); }

#define PTOK if (outfp) fprintf(outfp, "%s", yytext)
%}

%s GENERATOR
%s GENERATED

%%
<INITIAL>"/* <generator matexpr"           { PTOK; 
                                             BEGIN GENERATOR; 
                                             return GENBLOCK; }
<INITIAL>[ \t]*"/* <generated matexpr> */" { BEGIN GENERATED; }
<INITIAL>\n { PTOK; ++linenum; linepos = 0; }
<INITIAL>.  { PTOK; ++linepos; }


<GENERATOR>"*/"[ \t\r]*\n  { 
    PTOK;
    generate(outfp); 
    ++linenum; linepos = 0;
    print_cpp_line(linenum);
    BEGIN 0; 
}

<GENERATOR>"//"[^\n]*\n  { PTOK; ++linenum; }

<GENERATOR>"'"     { PTOK; return PRIME; }
<GENERATOR>"+="    { PTOK; return ADDTO; }

<GENERATOR>"complex"   { PTOK; return COMPLEX;   }
<GENERATOR>"symmetric" { PTOK; return SYMMETRIC; }
<GENERATOR>"input"     { PTOK; return INPUT;     }
<GENERATOR>"output"    { PTOK; return OUTPUT;    }
<GENERATOR>"inout"     { PTOK; return INOUT;     }
<GENERATOR>"function"  { PTOK; return FUNCTION;  }

<GENERATOR>[_a-zA-Z][_a-zA-Z0-9]* { 
    PTOK;
    yylval.s = strdup(yytext); 
    return ID; 
}

<GENERATOR>[0-9]+ { 
    PTOK;
    yylval.s = strdup(yytext);
    return INUMBER; 
}

<GENERATOR>(([0-9]+[.][0-9]*)|([0-9]*[.][0-9]+))([eE][+-]?[0-9]+)? { 
    PTOK;
    yylval.s = strdup(yytext);
    return FNUMBER; 
}

<GENERATOR>\"[^"\n]*\" { 
    PTOK; 
    yylval.s = strdup(yytext+1);
    yylval.s[strlen(yylval.s)-1] = 0;
    return STRING;
}

<GENERATOR>[ \t\r]+ { PTOK; }
<GENERATOR>\n       { PTOK; ++linenum; }
<GENERATOR>.        { PTOK; return yytext[0]; }


<GENERATED>"/* </generated matexpr> */"[ \t\r]*\n { 
    BEGIN 0; 
    ++linenum; 
    linepos = 0; 
}

<GENERATED>\n { ++ linenum; }
<GENERATED>. ;

%%
