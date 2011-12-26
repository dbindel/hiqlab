function impose_beam_x(x1,x2,y1,y2)

  local l  = x2-x1
  local ym = (y1+y2)/2;

  local function M1(x)  return 1+x*(  x*(-3+x* 2)),    x*(-6+x* 6);  end
  local function M2(x)  return   x*(1+x*(-2+x   )),  1+x*(-4+x* 3);  end
  local function M3(x)  return   x*(  x*( 3+x*-2)),    x*( 6+x*-6);  end
  local function M4(x)  return   x*(  x*(-1+x   )),    x*(-2+x* 3);  end

  local function N1(x,y) 
    if (meshbetween(x,x1,x2) and meshbetween(y,y1,y2)) then
      return (l-x)/l, 0; 
    end
  end

  local function N2(x,y)
    if (meshbetween(x,x1,x2) and meshbetween(y,y1,y2)) then
      local M, dM = M1(x/l)
      return -dM*(y-ym)/l, M
    end
  end

  local function N3(x,y)
    if (meshbetween(x,x1,x2) and meshbetween(y,y1,y2)) then
      local M, dM = M2(x/l)
      return -dM*(y-ym)/l, M
    end
  end

  local function N4(x,y) 
    if (meshbetween(x,x1,x2) and meshbetween(y,y1,y2)) then
      return x/l, 0; 
    end
  end

  local function N5(x,y)
    if (meshbetween(x,x1,x2) and meshbetween(y,y1,y2)) then
      local M, dM = M3(x/l)
      return -dM*(y-ym)/l, M
    end
  end

  local function N6(x,y)
    if (meshbetween(x,x1,x2) and meshbetween(y,y1,y2)) then
      local M, dM = M4(x/l)
      return -dM*(y-ym)/l, M
    end
  end

  local idx1 = mesh:add_global(N1)
  local idy1 = mesh:add_global(N2)
  local idt1 = mesh:add_global(N3)
  local idx2 = mesh:add_global(N4)
  local idy2 = mesh:add_global(N5)
  local idt2 = mesh:add_global(N6)

  if nil then
    clamp_boundary(
      function(x,y) 
        return (meshbetween(x,x1,x2) and meshbetween(y,y1,y2))
      end, 'ux', 'uy')
  end

  return {idx1, idy1, idt1, idx2, idy2, idt2}

end
