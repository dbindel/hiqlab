include ../../make.inc

NP=2
MPIRUN=mpirun -np $(NP)

test_petsc:
	$(MPIRUN) $(SRCDIR)/lua_hybrid/hiqlab petsc.lua

test_trilinos:
	$(MPIRUN) $(SRCDIR)/lua_hybrid/hiqlab trilinos.lua

test_serial:
	$(SRCDIR)/lua/hiqlab serial.lua

