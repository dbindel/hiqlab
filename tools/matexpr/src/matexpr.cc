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

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



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




/* Copy the first part of user declarations.  */
#line 1 "matexpr.y"

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



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 43 "matexpr.y"
typedef union YYSTYPE {
    struct ASTNode* ast;
    char* s;
} YYSTYPE;
/* Line 190 of yacc.c.  */
#line 155 "matexpr.cc"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 213 of yacc.c.  */
#line 167 "matexpr.cc"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

# ifndef YYFREE
#  define YYFREE free
# endif
# ifndef YYMALLOC
#  define YYMALLOC malloc
# endif

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   else
#    define YYSTACK_ALLOC alloca
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (defined (YYSTYPE_IS_TRIVIAL) && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short int yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short int) + sizeof (YYSTYPE))			\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined (__GNUC__) && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short int yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  17
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   126

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  30
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  21
/* YYNRULES -- Number of rules. */
#define YYNRULES  62
/* YYNRULES -- Number of states. */
#define YYNSTATES  112

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   271

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      25,    26,    19,    18,    27,    17,     2,    20,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    16,    24,
       2,    23,    22,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    28,     2,    29,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    21
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned char yyprhs[] =
{
       0,     0,     3,     6,     7,     8,    13,    18,    23,    32,
      36,    40,    41,    43,    45,    47,    50,    53,    57,    59,
      61,    65,    69,    71,    76,    83,    89,    99,   101,   103,
     107,   111,   115,   119,   123,   126,   128,   131,   135,   137,
     139,   141,   143,   148,   152,   155,   157,   158,   162,   165,
     168,   170,   174,   177,   180,   182,   183,   187,   190,   193,
     195,   196,   200
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      31,     0,    -1,    32,    31,    -1,    -1,    -1,    10,    33,
      34,    22,    -1,    12,    23,    40,    24,    -1,    12,     8,
      40,    24,    -1,    11,    12,    25,    49,    26,    23,    40,
      24,    -1,    35,    36,    24,    -1,     6,    23,    13,    -1,
      -1,     3,    -1,     4,    -1,     5,    -1,     6,     3,    -1,
       6,     5,    -1,    36,    27,    37,    -1,    37,    -1,    38,
      -1,    38,    23,    40,    -1,    38,     8,    40,    -1,    12,
      -1,    12,    25,    14,    26,    -1,    12,    25,    14,    27,
      14,    26,    -1,    12,     7,    25,    14,    26,    -1,    12,
      28,    39,    29,    25,    14,    27,    14,    26,    -1,    12,
      -1,    14,    -1,    40,    16,    40,    -1,    40,    18,    40,
      -1,    40,    17,    40,    -1,    40,    19,    40,    -1,    40,
      20,    40,    -1,    17,    40,    -1,    41,    -1,    41,     9,
      -1,    25,    40,    26,    -1,    12,    -1,    14,    -1,    15,
      -1,    42,    -1,    12,    25,    47,    26,    -1,    28,    43,
      29,    -1,    45,    44,    -1,    45,    -1,    -1,    24,    45,
      44,    -1,    24,    45,    -1,    40,    46,    -1,    40,    -1,
      27,    40,    46,    -1,    27,    40,    -1,    40,    48,    -1,
      40,    -1,    -1,    27,    40,    48,    -1,    27,    40,    -1,
      12,    50,    -1,    12,    -1,    -1,    27,    12,    50,    -1,
      27,    12,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned char yyrline[] =
{
       0,    62,    62,    62,    65,    65,    68,    74,    81,    86,
      89,    93,    96,    97,    98,    99,   100,   103,   104,   107,
     108,   114,   123,   133,   143,   154,   163,   177,   178,   181,
     182,   183,   184,   185,   186,   188,   191,   192,   193,   194,
     196,   197,   198,   204,   207,   208,   209,   212,   213,   216,
     217,   220,   221,   224,   225,   226,   229,   230,   233,   234,
     235,   238,   239
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "INPUT", "OUTPUT", "INOUT", "COMPLEX",
  "SYMMETRIC", "ADDTO", "PRIME", "GENBLOCK", "FUNCTION", "ID", "STRING",
  "INUMBER", "FNUMBER", "':'", "'-'", "'+'", "'*'", "'/'", "UMINUS", "'>'",
  "'='", "';'", "'('", "')'", "','", "'['", "']'", "$accept", "statements",
  "statement", "@1", "genparams", "iospec", "decllist", "decl", "declbase",
  "lda", "expr", "expr2", "matrix", "matdata", "rowrest", "rowdata",
  "colrest", "args", "argsrest", "formals", "formalsrest", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short int yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,    58,    45,    43,    42,
      47,   271,    62,    61,    59,    40,    41,    44,    91,    93
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    30,    31,    31,    33,    32,    32,    32,    32,    32,
      34,    34,    35,    35,    35,    35,    35,    36,    36,    37,
      37,    37,    38,    38,    38,    38,    38,    39,    39,    40,
      40,    40,    40,    40,    40,    40,    41,    41,    41,    41,
      41,    41,    41,    42,    43,    43,    43,    44,    44,    45,
      45,    46,    46,    47,    47,    47,    48,    48,    49,    49,
      49,    50,    50
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     2,     0,     0,     4,     4,     4,     8,     3,
       3,     0,     1,     1,     1,     2,     2,     3,     1,     1,
       3,     3,     1,     4,     6,     5,     9,     1,     1,     3,
       3,     3,     3,     3,     2,     1,     2,     3,     1,     1,
       1,     1,     4,     3,     2,     1,     0,     3,     2,     2,
       1,     3,     2,     2,     1,     0,     3,     2,     2,     1,
       0,     3,     2
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       3,    12,    13,    14,     0,     4,     0,     0,     0,     3,
       0,    15,    16,    11,     0,     0,     0,     1,     2,    22,
       0,    18,    19,     0,     0,    60,    38,    39,    40,     0,
       0,    46,     0,    35,    41,     0,     0,     0,     0,     9,
       0,     0,     0,     0,     5,    59,     0,    55,    34,     0,
      50,     0,    45,     0,     0,     0,     0,     0,     7,    36,
       6,     0,     0,    27,    28,     0,    17,    21,    20,    10,
       0,    58,     0,    54,     0,    37,     0,    49,    43,     0,
      44,    29,    31,    30,    32,    33,     0,    23,     0,     0,
      62,     0,     0,    53,    42,    52,    48,    25,     0,     0,
      61,     0,    57,    51,    47,    24,     0,     8,    56,     0,
       0,    26
};

/* YYDEFGOTO[NTERM-NUM]. */
static const yysigned_char yydefgoto[] =
{
      -1,     8,     9,    13,    24,    10,    20,    21,    22,    65,
      50,    33,    34,    51,    80,    52,    77,    74,    93,    46,
      71
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -16
static const yysigned_char yypact[] =
{
      59,   -16,   -16,   -16,    25,   -16,    -8,    -1,    31,    59,
       5,   -16,   -16,    19,     4,    -9,    -9,   -16,   -16,    -5,
      -6,   -16,     1,    20,    46,    60,    48,   -16,   -16,    -9,
      -9,    -9,    62,    65,   -16,    71,    50,    69,    43,   -16,
       5,    -9,    -9,    72,   -16,    57,    66,    -9,   -16,    34,
      17,    64,    70,    -9,    -9,    -9,    -9,    -9,   -16,   -16,
     -16,    87,    32,   -16,   -16,    73,   -16,    89,    89,   -16,
      91,   -16,    88,    29,    84,   -16,    -9,   -16,   -16,    -9,
     -16,    -7,    47,    47,   -16,   -16,    86,   -16,    99,    90,
      57,    -9,    -9,   -16,   -16,    17,    70,   -16,    92,   100,
     -16,    80,    29,   -16,   -16,   -16,    93,   -16,   -16,   102,
      95,   -16
};

/* YYPGOTO[NTERM-NUM].  */
static const yysigned_char yypgoto[] =
{
     -16,   108,   -16,   -16,   -16,   -16,   -16,    79,   -16,   -16,
     -15,   -16,   -16,   -16,    26,    44,    30,   -16,    22,   -16,
      36
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const unsigned char yytable[] =
{
      32,    35,    36,    26,    14,    27,    28,    15,    29,    41,
      54,    55,    56,    57,    48,    49,    30,    19,    39,    31,
      37,    40,    16,    38,    42,    23,    67,    68,    11,    25,
      12,    17,    73,    53,    54,    55,    56,    57,    81,    82,
      83,    84,    85,    43,    76,    53,    54,    55,    56,    57,
      53,    54,    55,    56,    57,    63,    92,    64,    87,    88,
      75,    95,     1,     2,     3,     4,    56,    57,    44,     5,
       6,     7,    45,    47,    59,    61,   101,   102,    53,    54,
      55,    56,    57,    62,    70,    69,    58,    53,    54,    55,
      56,    57,    72,    78,    79,    60,    53,    54,    55,    56,
      57,    86,    89,    90,   107,    53,    54,    55,    56,    57,
      94,    91,    97,    98,   106,    99,   110,    18,   105,    66,
     109,   111,   104,    96,   108,   103,   100
};

static const unsigned char yycheck[] =
{
      15,    16,     7,    12,    12,    14,    15,     8,    17,     8,
      17,    18,    19,    20,    29,    30,    25,    12,    24,    28,
      25,    27,    23,    28,    23,     6,    41,    42,     3,    25,
       5,     0,    47,    16,    17,    18,    19,    20,    53,    54,
      55,    56,    57,    23,    27,    16,    17,    18,    19,    20,
      16,    17,    18,    19,    20,    12,    27,    14,    26,    27,
      26,    76,     3,     4,     5,     6,    19,    20,    22,    10,
      11,    12,    12,    25,     9,    25,    91,    92,    16,    17,
      18,    19,    20,    14,    27,    13,    24,    16,    17,    18,
      19,    20,    26,    29,    24,    24,    16,    17,    18,    19,
      20,    14,    29,    12,    24,    16,    17,    18,    19,    20,
      26,    23,    26,    14,    14,    25,    14,     9,    26,    40,
      27,    26,    96,    79,   102,    95,    90
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,     3,     4,     5,     6,    10,    11,    12,    31,    32,
      35,     3,     5,    33,    12,     8,    23,     0,    31,    12,
      36,    37,    38,     6,    34,    25,    12,    14,    15,    17,
      25,    28,    40,    41,    42,    40,     7,    25,    28,    24,
      27,     8,    23,    23,    22,    12,    49,    25,    40,    40,
      40,    43,    45,    16,    17,    18,    19,    20,    24,     9,
      24,    25,    14,    12,    14,    39,    37,    40,    40,    13,
      27,    50,    26,    40,    47,    26,    27,    46,    29,    24,
      44,    40,    40,    40,    40,    40,    14,    26,    27,    29,
      12,    23,    27,    48,    26,    40,    45,    26,    14,    25,
      50,    40,    40,    46,    44,    26,    14,    24,    48,    27,
      14,    26
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");\
      YYERROR;							\
    }								\
while (0)


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (N)								\
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (0)
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
              (Loc).first_line, (Loc).first_column,	\
              (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr, 					\
                  Type, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short int *bottom, short int *top)
#else
static void
yy_stack_print (bottom, top)
    short int *bottom;
    short int *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);


# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
        break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short int yyssa[YYINITDEPTH];
  short int *yyss = yyssa;
  register short int *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;


  yyvsp[0] = yylval;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short int *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short int *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a look-ahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to look-ahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 4:
#line 65 "matexpr.y"
    {
      complexname = default_complexname;
    ;}
    break;

  case 6:
#line 68 "matexpr.y"
    { 
      ASTNode* stmt = new ASTNode(AST_ASSIGN, new ASTNode((yyvsp[-3].s)), (yyvsp[-1].ast)); 
      stmt->line = linenum;
      stmts.push_back(stmt);
      free((yyvsp[-3].s));
    ;}
    break;

  case 7:
#line 74 "matexpr.y"
    {
      ASTNode* stmt = new ASTNode(AST_ASSIGN, new ASTNode((yyvsp[-3].s)),
				  new ASTNode(AST_ADD, new ASTNode((yyvsp[-3].s)), (yyvsp[-1].ast)));
      stmt->line = linenum;
      stmts.push_back(stmt);
      free((yyvsp[-3].s));
    ;}
    break;

  case 8:
#line 81 "matexpr.y"
    {
      ASTNode* fn = new ASTNode(AST_FUNCTION, (yyvsp[-4].ast), (yyvsp[-1].ast));
      funtable[(yyvsp[-6].s)] = fn; // FIXME: Clean up old defs
      free((yyvsp[-6].s));
    ;}
    break;

  case 10:
#line 89 "matexpr.y"
    { 
      complexname = (yyvsp[0].s); 
      free((yyvsp[0].s)); 
    ;}
    break;

  case 12:
#line 96 "matexpr.y"
    { which_iospec = AST_INPUT;  ;}
    break;

  case 13:
#line 97 "matexpr.y"
    { which_iospec = AST_OUTPUT; ;}
    break;

  case 14:
#line 98 "matexpr.y"
    { which_iospec = AST_INOUT;  ;}
    break;

  case 15:
#line 99 "matexpr.y"
    { which_iospec = AST_INPUTZ;  ;}
    break;

  case 16:
#line 100 "matexpr.y"
    { which_iospec = AST_INOUTZ;  ;}
    break;

  case 19:
#line 107 "matexpr.y"
    { free((yyvsp[0].s)); ;}
    break;

  case 20:
#line 108 "matexpr.y"
    {
      ASTNode* stmt = new ASTNode(AST_ASSIGN, new ASTNode((yyvsp[-2].s)), (yyvsp[0].ast)); 
      stmt->line = linenum;
      stmts.push_back(stmt);
      free((yyvsp[-2].s));
    ;}
    break;

  case 21:
#line 114 "matexpr.y"
    {
      ASTNode* stmt = new ASTNode(AST_ASSIGN, new ASTNode((yyvsp[-2].s)), 
                                  new ASTNode(AST_ADD, new ASTNode((yyvsp[-2].s)), (yyvsp[0].ast))); 
      stmt->line = linenum;
      stmts.push_back(stmt);
      free((yyvsp[-2].s));
    ;}
    break;

  case 22:
#line 123 "matexpr.y"
    { 
      char ldabuf[16];
      ASTNode* stmt = new ASTNode(which_iospec, new ASTNode((yyvsp[0].s))); 
      stmt->line = linenum;
      stmt->m = 1;
      stmt->n = 1;
      stmt->lda = "1";
      stmt->array_flag = 0;
      stmts.push_back(stmt);
  ;}
    break;

  case 23:
#line 133 "matexpr.y"
    { 
      ASTNode* stmt = new ASTNode(which_iospec, new ASTNode((yyvsp[-3].s)));
      stmt->line = linenum;
      stmt->m = atoi((yyvsp[-1].s));
      stmt->n = 1;
      stmt->lda = (yyvsp[-1].s);
      stmt->array_flag = 1;
      stmts.push_back(stmt);
      free((yyvsp[-1].s));
  ;}
    break;

  case 24:
#line 143 "matexpr.y"
    {
      ASTNode* stmt = new ASTNode(which_iospec, new ASTNode((yyvsp[-5].s)));
      stmt->line = linenum;
      stmt->m = atoi((yyvsp[-3].s));
      stmt->n = atoi((yyvsp[-1].s));
      stmt->lda = (yyvsp[-3].s);
      stmt->array_flag = 1;
      stmts.push_back(stmt);
      free((yyvsp[-3].s));
      free((yyvsp[-1].s));
  ;}
    break;

  case 25:
#line 154 "matexpr.y"
    {
      ASTNode* stmt = new ASTNode(which_iospec, new ASTNode((yyvsp[-4].s)));
      stmt->line = linenum;
      stmt->m = stmt->n = atoi((yyvsp[-1].s));
      stmt->lda = (yyvsp[-1].s);
      stmt->array_flag = 2;
      stmts.push_back(stmt);
      free((yyvsp[-1].s));
  ;}
    break;

  case 26:
#line 163 "matexpr.y"
    {
      ASTNode* stmt = new ASTNode(which_iospec, new ASTNode((yyvsp[-8].s)));
      stmt->line = linenum;
      stmt->m = atoi((yyvsp[-3].s));
      stmt->n = atoi((yyvsp[-1].s));
      stmt->lda = (yyvsp[-6].s);
      stmt->array_flag = 1;
      stmts.push_back(stmt);
      free((yyvsp[-6].s));
      free((yyvsp[-3].s));
      free((yyvsp[-1].s));
  ;}
    break;

  case 27:
#line 177 "matexpr.y"
    { (yyval.s) = (yyvsp[0].s); ;}
    break;

  case 28:
#line 178 "matexpr.y"
    { (yyval.s) = (yyvsp[0].s); ;}
    break;

  case 29:
#line 181 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_RANGE, (yyvsp[-2].ast), (yyvsp[0].ast)); ;}
    break;

  case 30:
#line 182 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_ADD, (yyvsp[-2].ast), (yyvsp[0].ast)); ;}
    break;

  case 31:
#line 183 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_SUB, (yyvsp[-2].ast), (yyvsp[0].ast)); ;}
    break;

  case 32:
#line 184 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_MUL, (yyvsp[-2].ast), (yyvsp[0].ast)); ;}
    break;

  case 33:
#line 185 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_DIV, (yyvsp[-2].ast), (yyvsp[0].ast)); ;}
    break;

  case 34:
#line 187 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_NEG, (yyvsp[0].ast)); ;}
    break;

  case 36:
#line 191 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_TRANSP, (yyvsp[-1].ast)); ;}
    break;

  case 37:
#line 192 "matexpr.y"
    { (yyval.ast) = (yyvsp[-1].ast); ;}
    break;

  case 38:
#line 193 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_ID,    (yyvsp[0].s)); free((yyvsp[0].s)); ;}
    break;

  case 39:
#line 194 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_VALUE, (yyvsp[0].s)); free((yyvsp[0].s)); 
                    (yyval.ast)->name += ".0"; ;}
    break;

  case 40:
#line 196 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_VALUE, (yyvsp[0].s)); free((yyvsp[0].s)); ;}
    break;

  case 41:
#line 197 "matexpr.y"
    { (yyval.ast) = (yyvsp[0].ast); ;}
    break;

  case 42:
#line 198 "matexpr.y"
    { 
      (yyval.ast) = new ASTNode((yyvsp[-3].s), (yyvsp[-1].ast));
      free((yyvsp[-3].s));
  ;}
    break;

  case 43:
#line 204 "matexpr.y"
    { (yyval.ast) = (yyvsp[-1].ast); ;}
    break;

  case 44:
#line 207 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_VCAT, (yyvsp[-1].ast), (yyvsp[0].ast)); ;}
    break;

  case 45:
#line 208 "matexpr.y"
    { (yyval.ast) = (yyvsp[0].ast);   ;}
    break;

  case 46:
#line 209 "matexpr.y"
    { (yyval.ast) = NULL; ;}
    break;

  case 47:
#line 212 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_VCAT, (yyvsp[-1].ast), (yyvsp[0].ast)); ;}
    break;

  case 48:
#line 213 "matexpr.y"
    { (yyval.ast) = (yyvsp[0].ast); ;}
    break;

  case 49:
#line 216 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_HCAT, (yyvsp[-1].ast), (yyvsp[0].ast)); ;}
    break;

  case 51:
#line 220 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_HCAT, (yyvsp[-1].ast), (yyvsp[0].ast)); ;}
    break;

  case 52:
#line 221 "matexpr.y"
    { (yyval.ast) = (yyvsp[0].ast); ;}
    break;

  case 53:
#line 224 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_ARG, (yyvsp[-1].ast), (yyvsp[0].ast)); ;}
    break;

  case 54:
#line 225 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_ARG, (yyvsp[0].ast));     ;}
    break;

  case 55:
#line 226 "matexpr.y"
    { (yyval.ast) = NULL; ;}
    break;

  case 56:
#line 229 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_ARG, (yyvsp[-1].ast), (yyvsp[0].ast));   ;}
    break;

  case 57:
#line 230 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_ARG, (yyvsp[0].ast));       ;}
    break;

  case 58:
#line 233 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_ARG, new ASTNode((yyvsp[-1].s)), (yyvsp[0].ast)); ;}
    break;

  case 59:
#line 234 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_ARG, new ASTNode((yyvsp[0].s)));     ;}
    break;

  case 60:
#line 235 "matexpr.y"
    { (yyval.ast) = NULL; ;}
    break;

  case 61:
#line 238 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_ARG, new ASTNode((yyvsp[-1].s)), (yyvsp[0].ast)); ;}
    break;

  case 62:
#line 239 "matexpr.y"
    { (yyval.ast) = new ASTNode(AST_ARG, new ASTNode((yyvsp[0].s)));     ;}
    break;


    }

/* Line 1037 of yacc.c.  */
#line 1517 "matexpr.cc"

  yyvsp -= yylen;
  yyssp -= yylen;


  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (yychar);
	  const char* yyprefix;
	  char *yymsg;
	  int yyx;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 0;

	  yyprefix = ", expecting ";
	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		yysize += yystrlen (yyprefix) + yystrlen (yytname [yyx]);
		yycount += 1;
		if (yycount == 5)
		  {
		    yysize = 0;
		    break;
		  }
	      }
	  yysize += (sizeof ("syntax error, unexpected ")
		     + yystrlen (yytname[yytype]));
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yyprefix = ", expecting ";
		  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			yyp = yystpcpy (yyp, yyprefix);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yyprefix = " or ";
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("syntax error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* If at end of input, pop the error token,
	     then the rest of the stack, then return failure.  */
	  if (yychar == YYEOF)
	     for (;;)
	       {

		 YYPOPSTACK;
		 if (yyssp == yyss)
		   YYABORT;
		 yydestruct ("Error: popping",
                             yystos[*yyssp], yyvsp);
	       }
        }
      else
	{
	  yydestruct ("Error: discarding", yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

#ifdef __GNUC__
  /* Pacify GCC when the user code never invokes YYERROR and the label
     yyerrorlab therefore never appears in user code.  */
  if (0)
     goto yyerrorlab;
#endif

yyvsp -= yylen;
  yyssp -= yylen;
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping", yystos[yystate], yyvsp);
      YYPOPSTACK;
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token. */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yydestruct ("Error: discarding lookahead",
              yytoken, &yylval);
  yychar = YYEMPTY;
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 241 "matexpr.y"

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

