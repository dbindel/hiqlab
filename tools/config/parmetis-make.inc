
# Which compiler to use
CC = gcc


# What optimization level to use
OPTFLAGS = -O3 

# Include directories for the compiler
INCDIR = -I/opt/local/include

# What options to be used by the compiler
COPTIONS = 

# Which loader to use
LD = gcc

# In which directories to look for any additional libraries
LIBDIR = -L/opt/local/lib -L/usr/lib

# What additional libraries to link the programs with (eg., -lmpi)
#XTRALIBS = -lefence
#XTRALIBS = -ldmalloc

# What archiving to use
AR = /usr/bin/ar rv

# What to use for indexing the archive
#RANLIB = ranlib
RANLIB = ranlib

VERNUM = 
