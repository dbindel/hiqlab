mkoctfile --mex -g -I../corelib -I../numeric -I../element -I../io -I../lua -I/home/bindel/work/hiqlab/tools//UFconfig/ -I/home/bindel/work/hiqlab/tools//AMD/Include -I/home/bindel/work/hiqlab/tools//UMFPACK/Include -I/home/bindel/work/hiqlab/tools//lua/include -I/home/bindel/work/hiqlab/tools//tolua++/include meshmex.cc mexutil.cc meshstubs.cc luamatlab.cc luamex.cc /home/bindel/work/hiqlab/lib//libqlab.a -L/home/bindel/work/hiqlab/tools//tolua++/lib -ltolua -L/home/bindel/work/hiqlab/tools//lua/lib -llua -llualib /home/bindel/work/hiqlab/tools//ARPACK/libarpack.a /home/bindel/work/hiqlab/tools//UMFPACK/Lib/libumfpack.a /home/bindel/work/hiqlab/tools//AMD/Lib/libamd.a /usr/local/lib/liblapack.a -L/usr/local/lib -lgoto -lpthread -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lgfortranbegin -lgfortran -lm -lgcc_s -lm