all: otest ltest

otest:
	octave test_all.m

ltest:
	../bin/hiqlab test_qarray.lua
	../bin/hiqlab test_leigs.lua

clean:
	rm test_luaio.out
