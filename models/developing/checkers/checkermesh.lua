-- HiQLab
-- Copyright (c): Regents of the University of California


require 'common.lua'


-- Set default values for any unset parameters

w  = w  or 40e-6   -- Checker width
bw = bw or 4e-6    -- Border width (overlap region or beam width)
bl = bl or -bw     -- Beam length (if positive)
m1 = m1 or 4       -- Mesh density in checker border
m2 = m2 or 15      -- Mesh density in checker core
m3 = m3 or 4       -- Mesh density along beam length
nx = nx or 3       -- Number of checkers (X dir)
ny = ny or 3       -- Number of checkers (Y dir)
dx = dx or 0       -- Drive square (X coord)
dy = dy or 2       -- Drive square (Y coord)
sx = sx or 2       -- Sense square (X coord)
sy = sy or 0       -- Sense square (Y coord)
fs = fs or 0       -- Index of first square (zero if SE corner filled, 1 ow)

material = material or 'silicon2'


-- Mesh block spacing

meshtol = (w/(m1+2*m2)) / 100  -- Tolerance for tying
s    = w+bl                    -- Spacing between square corners
xmax = nx*s-bl                 -- x coord of leftmost square edge
ymax = ny*s-bl                 -- y coord of topmost square edge


-- Initialize mesh

order = 1
mesh  = Mesh:new(2)
elt   = mesh:PMLElastic2d_planestress(material)


-- Create a checker

function coordrange(i)
  local outerl = i*s
  local outerr = outerl + w
  local innerl = outerl + bw
  local innerr = outerr - bw
  return {outerl, innerl, innerr, outerr}
end

function checker(i,j)
  local xlist = coordrange(i)
  local ylist = coordrange(j)
  local mlist = {m1, m2, m1}
  mesh:blocks2dn(xlist, mlist, ylist, mlist, elt, order, dense)
end


-- Create beams from the ne or nw corners of the i,j square

function beam_ne(i,j)
  local x1,y1 = i*s+w, j*s+w
  local x2,y2 = x1+bl, y1+bl
  mesh:add_block_shape(m3,m1,elt,order, {x1,y1, x1-bw,y1, x2,y2, x2,y2+bw})
  mesh:add_block_shape(m3,m1,elt,order, {x1,y1-bw, x1,y1, x2+bw,y2, x2,y2})
end

function beam_nw(i,j)
  local x1,y1 = i*s, j*s+w
  local x2,y2 = x1-bl, y1+bl
  mesh:add_block_shape(m3,m1,elt,order, {x2,y2, x2,y2+bw, x1,y1, x1+bw,y1})
  mesh:add_block_shape(m3,m1,elt,order, {x2-bw,y2, x2,y2, x1,y1-bw, x1,y1})
end

function beam_sw(i,j)  beam_ne(i-1,j-1)  end
function beam_se(i,j)  beam_nw(i+1,j-1)  end

function beams(i,j)
  --
  -- Always put beams at NW and NE corners; only fill SE and SW when
  -- there is nothing to the SE (or to the SW) that would contribute.
  --
  if bl > 0 then
    beam_nw(i,j)
    beam_ne(i,j)
    if j == 0 or i == nx-1 then beam_se(i,j) end
    if j == 0 or i == 0    then beam_sw(i,j) end
  end
end


-- Mesh consists of checkers tied together

for i = 0,nx-1 do
  for j = 0,ny-1 do
    if mod(i+j,2) == fs then
      checker(i,j)
      beams(i,j)
    end
  end
end
mesh:tie()


-- Define boundary conditions

bc_list = {}
function add_bc(bc_func)
  bc_list[table.getn(bc_list)+1] = bc_func
end
mesh:set_bc(bc_list)


-- Build up a list of boundary conditions

add_bc(drive_checker)
for i = -1,nx do
  for j = -1,ny do
    if mod(i+j+2,2) == fs and        -- If we're on a legit checker...
        ((i == -1 or i == nx) or     -- and it's outside the boundary...
         (j == -1 or j == ny)) then

      -- Define a "ghost" rectangle where everything should be anchored
      -- if beam-coupled,   anchor everything inside the ghost checker;
      -- if corner-coupled, anchor everything inside the "core"
      --
      local l, r  = i*s, i*s+w
      local lo,hi = j*s, j*s+w
      if bl <= 0 then
        l,r,lo,hi = l+bw,r-bw,lo+bw,hi-bw
      end

      -- Add a predicate to the list to anchor inside the ghost
      --
      add_bc( 
        function(x,y)
          if meshbetween(x,l,r) and meshbetween(y,lo,hi) then return 'uu' end
        end)

    end
  end
end

function drive_checker(x,y)
  --
  -- Define forcing boundary conditions
  --
  local l,  r   = dx*s, dx*s+w
  local lo, hi  = dy*s, dy*s+w
  local el, er  = dx*s+w/3, dx*s+2*w/3
  local elo,ehi = dy*s+w/3, dy*s+2*w/3

  if     (mesheq(x,l)  and meshbetween(y,elo,ehi)) then return 'f ', -1
  elseif (mesheq(x,r)  and meshbetween(y,elo,ehi)) then return 'f ',  1
  elseif (mesheq(y,lo) and meshbetween(x,el, er )) then return ' f',  1
  elseif (mesheq(y,hi) and meshbetween(x,el, er )) then return ' f', -1
  end
end

function sense_checker(x,y)
  --
  -- Define sense pattern
  --
  local l,  r   = sx*s, sx*s+w
  local lo, hi  = sy*s, sy*s+w
  local el, er  = sx*s+w/3, sx*s+2*w/3
  local elo,ehi = sy*s+w/3, sy*s+2*w/3

  if     (mesheq(x,l)  and meshbetween(y,elo,ehi)) then return 'u ', -1
  elseif (mesheq(x,r)  and meshbetween(y,elo,ehi)) then return 'u ',  1
  elseif (mesheq(y,lo) and meshbetween(x,el, er )) then return ' u',  1
  elseif (mesheq(y,hi) and meshbetween(x,el, er )) then return ' u', -1
  end
end

