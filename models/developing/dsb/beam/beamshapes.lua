function impose_beam_x(x1,x2,y1,y2)

  local l  = x2-x1
  local ym = (y1+y2)/2;

  local function N1(x,y) return (l-x)/l, 0; end

  local function N2(x,y)
    x = x/l
    local M  = 1+x*(  x*(-3+x* 2))
    local dM = x*(-6+x* 6)
    return -dM*(y-ym)/l, M
  end

  local function N3(x,y)
    x = x/l
    local M  = x*(1+x*(-2+x   ))
    local dM = 1+x*(-4+x* 3);
    return -dM*(y-ym)/l, M
  end

  local function N4(x,y) return x/l, 0; end

  local function N5(x,y)
    x = x/l
    local M  = x*(  x*( 3+x*-2))
    local dM = x*( 6+x*-6)
    return -dM*(y-ym)/l, M
  end

  local function N6(x,y)
    x = x/l
    local M  = x*(  x*(-1+x   ))
    local dM = x*(-2+x* 3)
    return -dM*(y-ym)/l, M
  end

  local idx1 = mesh:add_global(N1)
  local idy1 = mesh:add_global(N2)
  local idt1 = mesh:add_global(N3)
  local idx2 = mesh:add_global(N4)
  local idy2 = mesh:add_global(N5)
  local idt2 = mesh:add_global(N6)

  clamp_boundary(
    function(x,y) 
      return (meshgeq(x,x1,1e-8) and meshleq(x,x2,1e-8) and
              meshgeq(y,y1,1e-8) and meshleq(y,y2,1e-8))
    end, 'ux', 'uy')

  return {idx1, idy1, idt1, idx2, idy2, idt2}

end
