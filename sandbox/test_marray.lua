-- Test Lua QArray interface

-- HiQLab
-- Copyright (c): Regents of the University of California
-- $Id: test_marray.lua,v 1.2 2004/10/18 21:07:16 dbindel Exp $

-- A = CSCMatrix
-- F = QArray

m,n = F:size()
x = QArray:new(m,n,2)
LU = UMFMatrix:new(A);
LU:solve(x,F)
LU:delete();
