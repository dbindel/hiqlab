/*
** Lua stand-alone interpreter
** See Copyright Notice in lua.h
**
** Modified by dbindel for HiQLab front-end
*/

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define lua_c

extern "C" {
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}

#include <tolua++.h>
#include "qhelperslua.h"
#include "qarraylua.h"
#include "qiarraylua.h"
#include "qlapacklua.h"
#include "cscmatrixlua.h"
#include "meshlua.h"
#include "elementlua.h"
#include "dxfilelua.h"
#include "shapeslua.h"
#include "leigslua.h"
#include "pmlmodelua.h"
#include "tedlinearlua.h"
#include "pzlinearlua.h"
#include "material_modellua.h"
#include "mesh_partitionlua.h"
#include "mesh_partitionerlua.h"
#include "mesh_add_blocklua.h"
#include "mesh_managerlua.h"

#include "qarray.h"
#include "qiarray.h"
#include "qlapack.h"
#include "mesh.h"
#include "dirstuff.h"
#include "hiqlab.h"

int    saved_argc;
char** saved_argv;

/*
** generic extra include file
*/
#ifdef LUA_USERCONFIG
#include LUA_USERCONFIG
#endif
#include "../../config.h"

#ifdef HAVE_FLTK
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>
#endif

#include "saconfig.cc"

// -- Petsc related
#include "petscksp.h"
#include "qpetsclua.h"
#include "qpassembly_petsclua.h"
#include "qpetsc_pclua.h"
#include "qpetsc_meshlua.h"
#include "qpetsc_jdqzlua.h"
#ifdef HAVE_HIQLAB_SLEPC
#include "qslepclua.h"
#include "slepceps.h"
#endif

/*
** definition of `isatty'
*/
#ifdef _POSIX_C_SOURCE
#include <unistd.h>
#define stdin_is_tty()        isatty(0)
#else
#define stdin_is_tty()        1  /* assume stdin is a tty */
#endif



#ifndef PROMPT
#define PROMPT                "> "
#endif


#ifndef PROMPT2
#define PROMPT2                ">> "
#endif

#ifndef PROGNAME
#define PROGNAME        "lua"
#endif

#ifndef lua_userinit
#define lua_userinit(L)                openstdlibs(L)
#endif


#ifndef LUA_EXTRALIBS
#define LUA_EXTRALIBS        /* empty */
#endif


static lua_State *L = NULL;

static const char *progname = PROGNAME;

static int hiqlab_defs(lua_State* state)
{
    tolua_open(state);
    tolua_qhelpers_open(state);
    tolua_qarray_open(state);
    tolua_qiarray_open(state);
    tolua_qlapack_open(state);
    tolua_cscmatrix_open(state);
    tolua_mesh_open(state);
    tolua_element_open(state);
    tolua_dxfile_open(state);
    tolua_shapes_open(state);
    tolua_leigs_open(state);
    tolua_pmlmode_open(state);
    tolua_tedlinear_open(state);
    tolua_pzlinear_open(state);
    tolua_material_model_open(state);
    tolua_mesh_partition_open(state);
    tolua_mesh_partitioner_open(state);
    tolua_mesh_add_block_open(state);
    tolua_mesh_manager_open(state);
    tolua_qpetsc_open(state);
    tolua_qpassembly_petsc_open(state);
    tolua_qpetsc_pc_open(state);
    tolua_qpetsc_mesh_open(state);
    tolua_qpetsc_jdqz_open(state);
#if HAVE_HIQLAB_SLEPC
    tolua_qslepc_open(state);
#endif
    lua_dirstuff_open(state);
    return 1;
}

static const luaL_reg lualibs[] = {
  {"base", luaopen_base},
  {"table", luaopen_table},
  {"io", luaopen_io},
  {"string", luaopen_string},
  {"math", luaopen_math},
  {"debug", luaopen_debug},
  {"loadlib", luaopen_loadlib},
  /* add your libraries here */
  {"hiqlab", hiqlab_defs},
  LUA_EXTRALIBS
  {NULL, NULL}
};



static void lstop (lua_State *l, lua_Debug *ar) {
  (void)ar;  /* unused arg. */
  lua_sethook(l, NULL, 0, 0);
  luaL_error(l, "interrupted!");
}


static void laction (int i) {
  signal(i, SIG_DFL); /* if another SIGINT happens before lstop,
                              terminate process (default action) */
  lua_sethook(L, lstop, LUA_MASKCALL | LUA_MASKRET | LUA_MASKCOUNT, 1);
}


static void print_usage (void) {
  fprintf(stderr,
  "usage: %s [options] [script [args]].\n"
  "Available options are:\n"
  "  -        execute stdin as a file\n"
  "  -e stat  execute string `stat'\n"
  "  -i       enter interactive mode after executing `script'\n"
  "  -l name  load and run library `name'\n"
  "  -v       show version information\n"
  "  --       stop handling options\n" ,
  progname);
}


static void l_message (const char *pname, const char *msg) {
  if (pname) fprintf(stderr, "%s: ", pname);
  fprintf(stderr, "%s\n", msg);
}


static int report (int status) {
  const char *msg;
  if (status) {
    msg = lua_tostring(L, -1);
    if (msg == NULL) msg = "(error with no message)";
    l_message(progname, msg);
    lua_pop(L, 1);
  }
  return status;
}


static int lcall (int narg, int clear) {
  int status;
  int base = lua_gettop(L) - narg;  /* function index */
  lua_pushliteral(L, "_TRACEBACK");
  lua_rawget(L, LUA_GLOBALSINDEX);  /* get traceback function */
  lua_insert(L, base);  /* put it under chunk and args */
  signal(SIGINT, laction);
  status = lua_pcall(L, narg, (clear ? 0 : LUA_MULTRET), base);
  signal(SIGINT, SIG_DFL);
  lua_remove(L, base);  /* remove traceback function */
  return status;
}


static void print_version (void) {
  l_message(NULL, "-------------------------------------------------------");
  l_message(NULL, HIQ_MESSAGE);
  l_message(NULL, LUA_VERSION "  " LUA_COPYRIGHT);
  l_message(NULL, "-------------------------------------------------------\n");
}


static void getargs (char *argv[], int n) {
  int i;
  lua_newtable(L);
  for (i=0; argv[i]; i++) {
    lua_pushnumber(L, i - n);
    lua_pushstring(L, argv[i]);
    lua_rawset(L, -3);
  }
  /* arg.n = maximum index in table `arg' */
  lua_pushliteral(L, "n");
  lua_pushnumber(L, i-n-1);
  lua_rawset(L, -3);
}


static int docall (int status) {
  if (status == 0) status = lcall(0, 1);
  return report(status);
}


static int file_input (const char *name) {
  return docall(luaL_loadfile(L, name));
}


static int dostring (const char *s, const char *name) {
  return docall(luaL_loadbuffer(L, s, strlen(s), name));
}


static int load_file (const char *name) {
  lua_pushliteral(L, "require");
  lua_rawget(L, LUA_GLOBALSINDEX);
  if (!lua_isfunction(L, -1)) {  /* no `require' defined? */
    lua_pop(L, 1);
    return file_input(name);
  }
  else {
    lua_pushstring(L, name);
    return report(lcall(1, 1));
  }
}


/*
** these macros can be used to perform initialization and finalization
** for lua_saveline and lua_readline
*/
#ifndef lua_initline
#define lua_initline(L,pname)        /* empty */
#endif

#ifndef lua_exitline
#define lua_exitline(L)                /* empty */
#endif


/*
** this macro can be used by some `history' system to save lines
** read in manual input
*/
#ifndef lua_saveline
#define lua_saveline(L,line)        /* empty */
#endif


/*
** this macro defines a function to show the prompt and reads the
** next line for manual input
*/
#ifndef lua_readline
#define lua_readline(L,prompt)                readline(L,prompt)

/* maximum length of an input line */
#ifndef MAXINPUT
#define MAXINPUT        512
#endif


static int readline (lua_State *l, const char *prompt) {
  static char buffer[MAXINPUT];
  if (prompt) {
    fputs(prompt, stdout);
    fflush(stdout);
  }
  if (fgets(buffer, sizeof(buffer), stdin) == NULL)
    return 0;  /* read fails */
  else {
    lua_pushstring(l, buffer);
    return 1;
  }
}

#endif


static const char *get_prompt (int firstline) {
  const char *p = NULL;
  lua_pushstring(L, firstline ? "_PROMPT" : "_PROMPT2");
  lua_rawget(L, LUA_GLOBALSINDEX);
  p = lua_tostring(L, -1);
  if (p == NULL) p = (firstline ? PROMPT : PROMPT2);
  lua_pop(L, 1);  /* remove global */
  return p;
}


static int incomplete (int status) {
  if (status == LUA_ERRSYNTAX &&
         strstr(lua_tostring(L, -1), "near `<eof>'") != NULL) {
    lua_pop(L, 1);
    return 1;
  }
  else
    return 0;
}


static int load_string (void) {
  int status;
  lua_settop(L, 0);
  if (lua_readline(L, get_prompt(1)) == 0)  /* no input? */
    return -1;
  if (lua_tostring(L, -1)[0] == '=') {  /* line starts with `=' ? */
    lua_pushfstring(L, "return %s", lua_tostring(L, -1)+1);/* `=' -> `return' */
    lua_remove(L, -2);  /* remove original line */
  }
  if (lua_tostring(L, -1)[0] == '!') {  /* line starts with '!' ? */
    system(lua_tostring(L,-1)+1);    /* Execute */
    status = luaL_loadbuffer(L, "", 0, "=stdin");
  } else {
    for (;;) {  /* repeat until gets a complete line */
      status = luaL_loadbuffer(L, lua_tostring(L, 1), lua_strlen(L, 1), "=stdin");
      if (!incomplete(status)) break;  /* cannot try to add lines? */
      if (lua_readline(L, get_prompt(0)) == 0)  /* no more input? */
        return -1;
      lua_concat(L, lua_gettop(L));  /* join lines */
    }
  }
  lua_saveline(L, lua_tostring(L, 1));
  lua_remove(L, 1);  /* remove line */
  return status;
}


static void manual_input (void) {
  int status;
  const char *oldprogname = progname;
  progname = NULL;
  lua_initline(L, PROGNAME); /* progname may contain a path, so use PROGNAME */
  while ((status = load_string()) != -1) {
    if (status == 0) status = lcall(0, 0);
    report(status);
    if (status == 0 && lua_gettop(L) > 0) {  /* any result to print? */
      lua_getglobal(L, "print");
      lua_insert(L, 1);
      if (lua_pcall(L, lua_gettop(L)-1, 0, 0) != 0)
        l_message(progname, lua_pushfstring(L, "error calling `print' (%s)",
                                               lua_tostring(L, -1)));
    }
  }
  lua_settop(L, 0);  /* clear stack */
  fputs("\n", stdout);
  lua_exitline(L);
  progname = oldprogname;
}


static void openstdlibs (lua_State *l) {
  const luaL_reg *lib = lualibs;
  for (; lib->func; lib++) {
    lib->func(l);  /* open library */
    lua_settop(l, 0);  /* discard any results */
  }
}


static int handle_luainit (void) {
  const char *init = getenv("HIQ_INIT");
  if (init == NULL)
    return 0;  /* status OK */
  else
    return file_input(init);
}


struct Smain {
    char filename[512];
    int interactive;
    int status;
};


static int pmain (lua_State *l) {
  struct Smain *s = (struct Smain *)lua_touserdata(l, 1);
  int status;
  int interactive = 0;
  L = l;
  lua_userinit(l);  /* open libraries */
  status = handle_luainit();
  if (s->filename[0])
      status = file_input(s->filename);
  if (s->interactive)
      manual_input();
  s->status = status;
  return 0;
}


int main(int argc, char *argv[])
{
    PetscErrorCode ierr;
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
#if HAVE_HIQLAB_SLEPC
    SlepcInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
#endif
    int gdbstate=0;
    while(gdbstate){}

    int status;
    struct Smain s;
    lua_State *l = lua_open();  // create state
    if (l == NULL) {
      l_message(argv[0], "cannot create state: not enough memory");
      return EXIT_FAILURE;
    }

    s.interactive = 0;
    s.filename[0] = '\0';
    ierr = PetscOptionsGetInt(PETSC_NULL, "-i", &(s.interactive), PETSC_NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetString(PETSC_NULL, "-f", s.filename, sizeof(s.filename), PETSC_NULL); CHKERRQ(ierr);

    status = lua_cpcall(l, &pmain, &s);
    report(status);

    lua_close(l);
#if HAVE_HIQLAB_SLEPC
    ierr = SlepcFinalize();CHKERRQ(ierr);
#endif
    ierr = PetscFinalize();CHKERRQ(ierr);
    return (status || s.status) ? EXIT_FAILURE : EXIT_SUCCESS;
}
