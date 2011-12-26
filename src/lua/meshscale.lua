-- Define non-dimensionalizing versions of bound functions

dim_scales = {}


function get_dim_scale(ftype)
  return dim_scales[ftype] or 1
end


function Mesh:get_scale(name)
  if type(name) == "number" then
    return name
  elseif dim_scales then
    return dim_scales[name] or 1
  else
    return 1
  end
end


function Mesh:set_scaling_vector()
  self:clear_scaling_vector()

  if dim_scales.vars then
    for i,s in dim_scales.vars do
      if mod(i,2) == 1 then       
        self:set_nodal_u_scale((i-1)/2, self:get_scale(s))
      else
        self:set_nodal_f_scale((i-2)/2, self:get_scale(s))
      end
    end
  end

  if dim_scales.elements then
    for i,elts in dim_scales.elements do
      for j,s in elts do
        if mod(j,2) == 1 then
          self:set_d1(self:ibranch((j-1)/2,i-1), self:get_scale(s))
        else
          self:set_d2(self:ibranch((j-2)/2,i-1), self:get_scale(s))
        end
      end
    end
  end

  if dim_scales.shapeg then
    for i,s in dim_scales.shapeg do
      if mod(i,2) == 1 then
        self:set_d1(self:iglobal((i-1)/2), self:get_scale(s))
      else
        self:set_d2(self:iglobal((i-2)/2), self:get_scale(s))
      end
    end
  end

end


local add_element2 = Mesh.add_element
function Mesh:add_element(e, etype, nen, n)
  n   = n or 1
  nen = nen or table.getn(e)
  local eltnum = add_element2(self, e, etype, nen, n) 
  if type(etype.aux_scales)=='table' then
    Mesh.set_aux_slots(self,eltnum,unpack(etype.aux_scales))
  end
  return eltnum
end


function Mesh:set_aux_slots(eltnum,...)

  -- Ensure that the dim_scales.elements array is set up
  dim_scales          = dim_scales or {}
  dim_scales.elements = dim_scales.elements or {}
  
  local aux_table = {}
  for i = 1,arg.n,2 do
    aux_table[i  ] = arg[i  ]
    aux_table[i+1] = arg[i+1]
  end
  dim_scales.elements[eltnum+1] = aux_table
end


-- Define and redefine functions for global shape functions
local add_global1 = Mesh.add_global
function Mesh:add_global(func, s1, s2)

  -- Ensure that the global_funcs array is set up
  -- Ensure that the dim_scales.shapeg array is set up
  self.global_funcs = self.global_funcs or {}
  dim_scales        = dim_scales or {}
  dim_scales.shapeg = dim_scales.shapeg or {}

  -- Lua based tables should be 1-based 
  local idg = add_global1(self, 1)
  self.global_funcs[idg+1]     = func

  if s1 and s2 then
    dim_scales.shapeg[2*idg+1] = s1
    dim_scales.shapeg[2*idg+2] = s2
  end
  return idg
end 


local add_global2 = Mesh.add_global
function Mesh:add_global(func, ...)
  local idg
  local numg = 0
  if type(func)=='table' then
    numg = table.getn(func)
    if (arg.n > 0) and (numg ~= 2*arg.n) then
      error('Number of functions and nondim variables must match')
    end
    for i = 1,numg do
      idg = add_global2(self, func[i],arg[2*i-1],arg[2*i])
    end
  else
    numg = 1
    idg  = add_global2(self, func,arg[1],arg[2])
  end
  idg = idg - numg + 1
  return idg
end

