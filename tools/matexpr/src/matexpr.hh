/* A Bison parser, made by GNU Bison 2.0.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     INPUT = 258,
     OUTPUT = 259,
     INOUT = 260,
     COMPLEX = 261,
     SYMMETRIC = 262,
     ADDTO = 263,
     PRIME = 264,
     GENBLOCK = 265,
     FUNCTION = 266,
     ID = 267,
     STRING = 268,
     INUMBER = 269,
     FNUMBER = 270,
     UMINUS = 271
   };
#endif
#define INPUT 258
#define OUTPUT 259
#define INOUT 260
#define COMPLEX 261
#define SYMMETRIC 262
#define ADDTO 263
#define PRIME 264
#define GENBLOCK 265
#define FUNCTION 266
#define ID 267
#define STRING 268
#define INUMBER 269
#define FNUMBER 270
#define UMINUS 271




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 43 "matexpr.y"
typedef union YYSTYPE {
    struct ASTNode* ast;
    char* s;
} YYSTYPE;
/* Line 1318 of yacc.c.  */
#line 74 "matexpr.hh"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



