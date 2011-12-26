-- HiQLab
-- Copyright (c): Regents of the University of California
-- $Id: beammesh2.lua,v 1.2 2006/06/19 16:59:48 dbindel Exp $

-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh  = Mesh:new(2)

-- Define approx size of elements and order
dense = 1         -- Approximate element size (for block generator)
order = 3         -- Order of elements

-- Define size of block
l = 10
w = 10

-- Define element type
mat   = mesh:PMLElastic2d(1, nu, 1, 0);

-- Define mesh using block command
mesh:blocks2d( { 0,l }, { 0,w }, mat )

