
--
-- Helper functions -- set relevant fields in the mesh object
--   for BC functions or function tables
--
function Mesh:set_bc(f)           self.bcfunc = f         end
function Mesh:set_globals_bc(f)   self.shapegbc = f       end
function Mesh:set_elements_bc(f)  self.elementsbc = f     end
function Mesh:set_globals(f)      self.global_funcs = f   end


--
-- Iterate over the output of a BC function
-- 
function bc_iterator(s, ...)
  local i = 0
  local j = 0
  local space_char = string.byte(" ")
  local slength = s and string['len'](s)
  return function()
    while s and i < slength do
      i = i + 1
      if string.byte(s,i) ~= space_char then
        j = j + 1
        return i, arg[j] or 0, string.sub(s,i,i)
      end
    end
  end
end


--
-- Convert a method that takes one argument into a method that
-- can optionally take a table of arguments
--
local function table_wrap(fname)
  local func = Mesh[fname]
  Mesh[fname] = function(self, fs)
    if type(fs) == "table" then
      for i,f in fs do func(self, f) end
    elseif fs then
      func(self, fs)
    end
  end
end


--------------------------------------------------------------------


--
-- Apply nodal BCs specified by bcfunc
--
function Mesh:form_nodal_bcs(bcfunc)
  for j = 1,self:numnp() do
    for i,value,type in bc_iterator(bcfunc(self:x(j-1))) do
      self:set_bcode(i-1,j-1,type)
      self:set_bv(i-1,j-1,value)
    end
  end
end


--
-- Apply element dof BCs specified by bcfunc
--
function Mesh:form_element_bcs(bcfunc)
  for j = 1,self:numelt() do
    for i,value,type in bc_iterator(bcfunc(j-1)) do
      self:set_bcode(self:ibranch(i-1,j-1),type)
      self:set_bv(self:ibranch(i-1,j-1),value)
    end
  end
end


--
-- Apply global dof BCs specified by bcfunc
--
function Mesh:form_global_bcs(bcfunc)
  for j = 1,self:numglobals() do
    local type, value = bcfunc(j-1)
    if type then
      self:set_bcode(self:iglobal(j-1),type)
      self:set_bv(self:iglobal(j-1),value or 0)
    end
  end
end


--
-- Allow BC applicators to take a list of functions
--
table_wrap 'form_nodal_bcs'
table_wrap 'form_global_bcs'
table_wrap 'form_element_bcs'


--
-- Apply the nodal, element, and global BCs, then set zero BCs on dofs
-- covered by some global.  This function can be replaced by a user function,
-- if that's more convenient!
--
function Mesh:form_bcs()
  self:form_nodal_bcs(self.bcfunc)
  self:form_global_bcs(self.shapegbc)
  self:form_element_bcs(self.elementsbc)
  self:apply_bc_shapeg_const()
end


--------------------------------------------------------------------

--
-- Compute contribution of BC function sensefunc evaluated on nodes.
--
function make_nodal_func(sensefunc, vtype)
  return function(self,v)
    for j = 1,self:numnp() do
      for i,value,type in bc_iterator(sensefunc(self:x(j-1))) do
        if type == vtype then
          v:set(self:inode(i-1,j-1),value)
        end
      end
    end
  end
end


--
-- Compute contribution of BC function sensefunc evaluated on elements
--
function make_element_func(sensefunc, vtype)
  return function(self,v)
    for j = 1,self:numelt() do
      for i,value,type in bc_iterator(sensefunc(j-1)) do
        if type == vtype then
          v:set(self:ibranch(i-1,j-1),value)
        end
      end
    end
  end
end


--
-- Compute contribution of BC function sensefunc evaluated on globals
--
function make_global_func(sensefunc, vtype)
  return function(self,v)
    for j = 1,self:numglobals() do
      local type, value = sensefunc(j-1)
      if type == vtype then
        v:set(self:iglobal(j-1),value)
      end
    end
  end
end


--------------------------------------------------------------------


--
-- Evaluate shape functions j at node i
--
function Mesh:get_shapeg(i,j)
  local ndf   = self:get_ndf()
  local nodei = floor(i/ndf)
  local dofi  = mod(i,ndf)
  if nodei < self:numnp() and self.global_funcs and self.global_funcs[j+1] then
    local results = {self.global_funcs[j+1](self:x(nodei))}
    return results[dofi+1] or 0
  else
    return 0
  end
end


--
-- Evaluate c * shape function j and assemble into v
--
function Mesh:get_shapeg_vec(v, c, j)
  if self.global_funcs and self.global_funcs[j+1] then
    for i = 1,self:numnp() do
      local results = {self.global_funcs[j+1](self:x(i-1))}
      for l,val in results do
        v:add(self:inode(l-1,i-1),c*val)
      end
    end
  end
end


--
-- Set zero displacement BCs for dofs in the support of some global shape
--
function Mesh:apply_bc_shapeg_const()
  if not self.global_funcs then return end
  for j,shape in self.global_funcs do
    for i=1,self:numnp() do
      local results = {shape(self:x(i-1))}
      for l,val in results do
        if val ~= 0 then
          self:set_bcode(l-1,i-1, 'u')
          self:set_bv(l-1,i-1, 0)
        end
      end
    end
  end
end
