-- HiQLab
-- Copyright (c): Regents of the University of California
-- $Id: beammesh.lua,v 1.2 2006/05/02 04:16:51 tkoyama Exp $

require 'common.lua'

l = 10e-6         -- Beam length
w = 2e-6          -- Beam width
dense = 0.5e-6    -- Approximate element size (for block generator)
-- dense = 1e-6    -- Approximate element size (for block generator)
order = 2         -- Order of elements

mesh  = Mesh:new(2)
mat   = make_material_e(mesh, 'silicon2', 'planestrain')
mesh:blocks2d( { 0, l }, { -w/2.0, w/2.0 }, mat )

mesh:set_bc(function(x,y)
  if x == 0 then return 'uu', 0, 0; end
end)

function force_tip(x,y)
  if x == l and y == 0 then return ' f', 1; end
end
