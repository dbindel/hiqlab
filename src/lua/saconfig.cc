/* saconfig.c -- configuration for stand-alone Lua interpreter
*
* #define LUA_USERCONFIG to this file
*
* Here are the features that can be customized using #define:
*
*** Line editing and history:
*   #define USE_READLINE to use the GNU readline library.
*
*   To use another library for this, use the code below as a start.
*   Make sure you #define lua_{read,save,init,exit}line accordingly.
*   If you do not #define lua_readline, you'll get a version based on fgets
*   that uses a static buffer of size MAXINPUT.
*
*
*** Static Lua libraries to be loaded at startup:
*   #define lua_userinit(L) to a Lua function that loads libraries; typically
*        #define lua_userinit(L)        openstdlibs(L);myinit(L)
*   or
*        #define lua_userinit(L)        myinit(L)
*
*   Another way is to add the prototypes of the init functions here and
*   #define LUA_EXTRALIBS accordingly. For example,
*        #define LUA_EXTRALIBS {"mylib","luaopen_mylib"},
*   Note the ending comma!
*
*
*** Prompts:
*   The stand-alone Lua interpreter uses two prompts: PROMPT and PROMPT2.
*   PROMPT is the primary prompt, shown when the intepreter is ready to receive
*   a new statement. PROMPT2 is the secondary prompt, shown while a statement
*   is being entered but is still incomplete.
*
*
*** Program name:
*   Error messages usually show argv[0] as a program name. In systems that do
*   not give a valid string as argv[0], error messages show PROGNAME instead.
*
*
*/

#ifdef HAVE_READLINE
/*
* This section implements lua_xxxxline for lua.c using the GNU readline
* and history libraries or compatible replacements.
*
* It has been successfully tested with:
*
* GNU    readline 2.2.1  (1998-07-17)
* GNU    readline 4.0    (1999-02-18) [harmless compiler warning]
* GNU    readline 4.3    (2002-07-16)
* NETBSD libedit  2.6.5  (2002-03-25)
* NETBSD libedit  2.6.9  (2004-05-01)
*/

#define lua_initline        myinitline
#define lua_exitline        myexitline
#define lua_readline        myreadline
#define lua_saveline        mysaveline

#include <ctype.h>

#ifdef HAVE_READLINE_READLINE_H
#  include <readline/readline.h>
#else
#  ifdef HAVE_READLINE_H
#    include <readline.h>
#  endif
#endif

#ifdef HAVE_READLINE_HISTORY_H
#  include <readline/history.h>
#endif

/* Environment variable names for the history file and the history size */
#ifndef LUA_HISTORY_ENV
#define LUA_HISTORY_ENV                "HIQ_HISTORY"
#endif

#ifndef LUA_HISTSIZE_ENV
#define LUA_HISTSIZE_ENV        "HIQ_HISTSIZE"
#endif

static char *myhist;
static int myhistsize;

static lua_State *myL;        /* readline does not pass user data to callbacks */

/* Read a line from the terminal with line editing */
static int myreadline(lua_State *L, const char *prompt)
{
  char *s;
  if (!(s = readline(prompt))) return 0;
  lua_pushstring(L, s);
  lua_pushliteral(L, "\n");
  lua_concat(L, 2);
  free(s);
  return 1;
}

/* Add a line to the history */
static void mysaveline(lua_State *L, const char *s)
{
#ifdef HAVE_READLINE_HISTORY
  const char *p;
  for (p = s; isspace(*p); p++) ;
  if (*p) {
    size_t n = strlen(s)-1;
    if (s[n] != '\n') {
      add_history(s);
    } else {
      lua_pushlstring(L, s, n);
      s = lua_tostring(L, -1);
      add_history(s);
      lua_pop(L, 1);
    }
  }
#endif
}

/* Reserved lua keywords */
static const char * const reskeywords[] = {
  "and", "break", "do", "else", "elseif", "end", "false",
  "for", "function", "if", "in", "local", "nil", "not", "or",
  "repeat", "return", "then", "true", "until", "while", NULL
};

static int valididentifier(const char *s)
{
  if (!(isalpha(*s) || *s == '_')) return 0;
  for (s++; *s; s++) if (!(isalpha(*s) || isdigit(*s) || *s == '_')) return 0;
  return 1;
}

/* Dynamically resizable match list */
typedef struct {
  char **list;
  size_t idx, allocated, matchlen;
} dmlist;

/* Add prefix + string + suffix to list and compute common prefix */
static int dmadd(dmlist *ml, const char *p, size_t pn, const char *s, int suf)
{
  char *t = NULL;

  if (ml->idx+1 >= ml->allocated &&
      !(ml->list = (char**) realloc(ml->list, sizeof(char *)*(ml->allocated += 32))))
    return -1;

  if (s) {
    size_t n = strlen(s);
    if (!(t = (char *)malloc(sizeof(char)*(pn+n+(suf?2:1))))) return 1;
    memcpy(t, p, pn);
    memcpy(t+pn, s, n);
    n += pn;
    t[n] = suf;
    if (suf) t[++n] = '\0';

    if (ml->idx == 0) {
      ml->matchlen = n;
    } else {
      size_t i;
      for (i = 0; i < ml->matchlen && i < n && ml->list[1][i] == t[i]; i++) ;
      ml->matchlen = i;        /* Set matchlen to common prefix */
    }
  }

  ml->list[++ml->idx] = t;
  return 0;
}

/* Get __index field of metatable of object on top of stack */
static int getmetaindex(lua_State *L)
{
  if (!lua_getmetatable(L, -1)) { lua_pop(L, 1); return 0; }
  lua_pushstring(L, "__index");
  lua_rawget(L, -2);
  lua_replace(L, -2);
  if (lua_isnil(L, -1) || lua_rawequal(L, -1, -2)) { lua_pop(L, 2); return 0; }
  lua_replace(L, -2);
  return 1;
} /* 1: obj -- val, 0: obj -- */

/* Get field from object on top of stack. Avoid calling metamethods */
static int safegetfield(lua_State *L, const char *s, size_t n)
{
  int i = 20; /* Avoid infinite metatable loops */
  do {
    if (lua_istable(L, -1)) {
      lua_pushlstring(L, s, n);
      lua_rawget(L, -2);
      if (!lua_isnil(L, -1)) { lua_replace(L, -2); return 1; }
      lua_pop(L, 1);
    }
  } while (--i > 0 && getmetaindex(L));
  lua_pop(L, 1);
  return 0;
} /* 1: obj -- val, 0: obj -- */

/* Completion function */
static char **mycomplete(const char *text, int start, int end)
{
  dmlist ml;
  const char *s;
  size_t i, n, dot;
  int savetop;

  if (!(text[0] == '\0' || isalpha(text[0]) || text[0] == '_')) return NULL;

  ml.list = NULL;
  ml.idx = ml.allocated = ml.matchlen = 0;

  savetop = lua_gettop(myL);
  lua_pushvalue(myL, LUA_GLOBALSINDEX);
  for (n = (size_t)(end-start), i = dot = 0; i < n; i++)
    if (text[i] == '.' || text[i] == ':') {
      if (!safegetfield(myL, text+dot, i-dot)) goto error; /* invalid prefix */
      dot = i+1; /* points to first char after dot/colon */
    }

  /* Add all matches against keywords if there is no dot/colon */
  if (dot == 0)
    for (i = 0; (s = reskeywords[i]) != NULL; i++)
      if (!strncmp(s, text, n) && dmadd(&ml, NULL, 0, s, ' ')) goto error;

  /* Add all valid matches from all tables/metatables */
  i = 20; /* Avoid infinite metatable loops */
  do {
    if (lua_istable(myL, -1))
      for (lua_pushnil(myL); lua_next(myL, -2); lua_pop(myL, 1))
        if (lua_type(myL, -2) == LUA_TSTRING) {
          s = lua_tostring(myL, -2);
          /* Only match names starting with '_' if explicitly requested */
          if (!strncmp(s, text+dot, n-dot) && valididentifier(s) &&
              (*s != '_' || text[dot] == '_')) {
            int suf = ' '; /* default suffix is a space */
            switch (lua_type(myL, -1)) {
            case LUA_TTABLE:        suf = '.'; break; /* No way to guess ':' */
            case LUA_TFUNCTION:        suf = '('; break;
            case LUA_TUSERDATA:
              if (lua_getmetatable(myL, -1)) { lua_pop(myL, 1); suf = ':'; }
              break;
            }
            if (dmadd(&ml, text, dot, s, suf)) goto error;
          }
        }
  } while (--i > 0 && getmetaindex(myL));

  if (ml.idx > 1) {
    /* list[0] holds the common prefix of all matches (may be "") */
    if (!(ml.list[0] = (char *)malloc(sizeof(char)*(ml.matchlen+1)))) {
error:
      lua_settop(myL, savetop);
      return NULL;
    }
    memcpy(ml.list[0], ml.list[1], ml.matchlen);
    ml.list[0][ml.matchlen] = '\0';
    /* Add the NULL list terminator */
    if (dmadd(&ml, NULL, 0, NULL, 0)) goto error;
  } else if (ml.idx == 1) {
    ml.list[0] = ml.list[1];                /* Only return common prefix */
    ml.list[1] = NULL;
  }

  lua_settop(myL, savetop);
  return ml.list;
}

/* Initialize library */
static void myinitline(lua_State *L, char *pname)
{
  char *s;

  myL = L;

  /* This allows for $if hiqlab ... $endif in ~/.inputrc */
  rl_readline_name = pname;
  /* Break words at every non-identifier character except '.' and ':' */
  rl_completer_word_break_characters =
    "\t\r\n !\"#$%&'()*+,-/;<=>?@[\\]^`{|}~";
  rl_completer_quote_characters = "\"'";
  rl_completion_append_character = '\0';
  rl_attempted_completion_function = mycomplete;
  rl_initialize();

#ifdef HAVE_READLINE_HISTORY
  /* Start using history, optionally set history size and load history file */
  using_history();
  if ((s = getenv(LUA_HISTSIZE_ENV)) &&
      (myhistsize = atoi(s))) stifle_history(myhistsize);
  if ((myhist = getenv(LUA_HISTORY_ENV))) read_history(myhist);
#endif
}

/* Finalize library */
static void myexitline(lua_State *L)
{
#ifdef HAVE_READLINE_HISTORY
  /* Optionally save history file */
  if (myhist) write_history(myhist);
#endif
}
#endif
