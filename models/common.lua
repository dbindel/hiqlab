-- HiQLab
-- Copyright (c): Regents of the University of California


-- Common function definitions

require 'materials.lua'
require 'components.lua'


-- Convert material name to table if needed
-- 
function get_material(mtype)
  if type(mtype) == 'string' then mtype = _G[mtype] end
  return mtype
end


-- Set characteristic scales for mechanical problem
-- 
function mech_nondim(mtype, cL)
  mtype = get_material(mtype);

  dim_scales    = dim_scales or {}
  dim_scales.M  = mtype.rho * cL^3
  dim_scales.L  = cL
  dim_scales.T  = cL/sqrt(mtype.E/mtype.rho)

  dim_scales.F  = dim_scales.M * dim_scales.L / dim_scales.T^2
  dim_scales.F2 = dim_scales.M                / dim_scales.T^2
  dim_scales.F1 = dim_scales.M / dim_scales.L / dim_scales.T^2
  return dim_scales
end


-- Set characteristic scales for TED problem
-- 
function ted_nondim(mtype, cL)
  mtype = get_material(mtype);
  local T0 = mtype.T0 or dim_scales.T0
  mech_nondim(mtype, cL)

  dim_scales.Th = T0 * mtype.at * mtype.E / mtype.rho / mtype.cp
  dim_scales.Qt = dim_scales.M * dim_scales.L^2/ dim_scales.T^3
  dim_scales.Qt2= dim_scales.M * dim_scales.L  / dim_scales.T^3
  dim_scales.Qt1= dim_scales.M                 / dim_scales.T^3
  return dim_scales
end


-- Set derived electrical scales
--
local function electrical_nondim(A)
  dim_scales.A  = A
  dim_scales.V  = dim_scales.M * dim_scales.L^2 / dim_scales.A / dim_scales.T^3
  dim_scales.Q  = dim_scales.A                  * dim_scales.T
  dim_scales.Q2 = dim_scales.A / dim_scales.L   * dim_scales.T
  dim_scales.Q1 = dim_scales.A / dim_scales.L^2 * dim_scales.T
  dim_scales.E  = dim_scales.V / dim_scales.L
  dim_scales.R  = dim_scales.V / dim_scales.A
  dim_scales.Li = dim_scales.V / dim_scales.A * dim_scales.T
  dim_scales.C  = dim_scales.Q / dim_scales.V
end


-- Set characteristic scales for piezo problem
-- 
function pz_nondim(mtype, cL)
  mtype = get_material(mtype);
  mech_nondim(mtype, cL)
  electrical_nondim(mtype.kds/mtype.d * cL^2 / dim_scales.T)
  return dim_scales
end


-- Set characteristic scales for electro-mech problem
-- 
function em_nondim(mtype, cL, eps)
  mtype = get_material(mtype);
  local eps = eps or mtype.eps or dim_scales.eps
  mech_nondim(mtype, cL)
  electrical_nondim(sqrt(eps*mtype.E*cL^3)/dim_scales.T)
  return dim_scales
end


---------------------------

-- Compute crystal axes for a cubic material given a wafer type and
-- the angle between the x axis and a projection of a material
-- axis onto the x-y plane
-- 
function wafer_orientation(wafer, angle)
  local axis1, axis2, axis

  -- Wafer orientation and crystal axes
  if wafer == '100' then       -- [100] wafer
    axis1 = { 1, 0, 0}
    axis2 = { 0, 1, 0}
  elseif wafer == '111' then   -- [111] wafer
    axis1 = {      sqrt(2.0/3.0),           0.0,   1.0/sqrt(3.0) }
    axis2 = {     -1.0/sqrt(6.0), 1.0/sqrt(2.0),   1.0/sqrt(3.0) }
  else
    error('Unknown wafer type: ' .. wafer)
  end

  -- Rotate crystal axes
  axis1 = { cos(angle)*axis1[1] - sin(angle)*axis1[2], 
            sin(angle)*axis1[1] + cos(angle)*axis1[2],
            axis1[3] }
  axis2 = { cos(angle)*axis2[1] - sin(angle)*axis2[2], 
            sin(angle)*axis2[1] + cos(angle)*axis2[2],
            axis2[3] }

  return axis1, axis2
end


-- Compute Lame parameters for an isotropic material
local function get_lame(mtype)
  local E      = mtype.E
  local nu     = mtype.nu
  local lambda = mtype.lambda or (E and nu and (nu*E/(1+nu)/(1-2*nu))) or
                 error('Could not compute modulus lambda')
  local mu     = mtype.mu     or (E and nu and (E/2/(1+nu))) or
                 error('Could not compute modulus nu')
  return lambda, mu
end


-- Set piezo materials
-- 
--  axis = [ axis1[1], axis1[2], axis1[3],
--           axis2[1], axis2[2], axis2[3] ]
--
local function wrap_pz(constructor, etype)
  local compute_D = {
    ['pstrain'] = piezo_elasticity_2D_strain,
    ['pstress'] = piezo_elasticity_2D_stress,
    ['2hd']     = piezo_elasticity_2HD_stress,
    ['3d']      = piezo_elasticity_3D
  }

  return function(self, mtype, wafer, angle)
    local axis1, axis2
    mtype = fill_piezo(get_material(mtype))
    if not angle then 
      local axis = wafer 
      axis1= { axis[1], axis[2], axis[3] }
      axis2= { axis[4], axis[5], axis[6] }
    else
      axis1,axis2 = wafer_orientation(wafer,angle) 
    end

    local Dmech = {}
    local pz = {}
    local kds = {}
    hex_elasticity_3D(Dmech,
                      { mtype.c11, mtype.c12, mtype.c13, mtype.c33, mtype.c55 },
                      axis1, axis2)
    hex_piezo_3D(pz, { mtype.d16, mtype.d31, mtype.d33 }, axis1, axis2)
    hex_dielectric_3D(kds, { mtype.kds1, mtype.kds3 }, axis1, axis2)

    local Db = {}
    compute_D[etype](Db, Dmech, pz, kds)
    return constructor(self, Db, mtype.rho or 0)
  end
end

Mesh.PMLElastic2d_pz_planestrain = wrap_pz(Mesh.PMLElastic2d_pz, 'pstrain')
Mesh.PMLElastic2d_pz_planestress = wrap_pz(Mesh.PMLElastic2d_pz, 'pstress')
Mesh.PMLElastic2hd_pz            = wrap_pz(Mesh.PMLElastic2d_pz, '2hd')
Mesh.PMLElastic3d_pz             = wrap_pz(Mesh.PMLElastic2d_pz, '3d')


-- Set TED materials
-- 
local function wrap_te(constructor, etype)
  local compute_D = {
    ['pstrain'] = thermoelasticity_2D_strain,
    ['pstress'] = thermoelasticity_2D_stress,
    ['axis']    = thermoelasticity_axis,
    ['3d']      = thermoelasticity_3D
  }

  local compute_D_cubic = {
    ['pstrain'] = cubic_thermoelasticity_2D_strain,
    ['pstress'] = cubic_thermoelasticity_2D_stress,
    ['3d']      = cubic_thermoelasticity_3D
  }

  return function(self, mtype, wafer, angle)
    mtype    = fill_mech(get_material(mtype))
    mtype.T0 = mtype.T0 or dim_scales.T0
    local Db = {}
    if wafer and angle then  
      local axis1, axis2 = wafer_orientation(wafer, angle)
      local f = compute_D_cubic[etype] or
                error('Cannot make te cubic ' .. etype .. ' material')
      f(Db, mtype.lambda, mtype.mu, mtype.alpha, axis1, axis2, mtype.at)
    else
      local lambda, mu = get_lame(mtype)
      compute_D[etype](Db, lambda, mu, mtype.at)
    end
    return constructor(self, Db, mtype.rho, mtype.at, mtype.cp, 
                       mtype.kt, mtype.T0)
  end
end

Mesh.PMLElastic2d_te_planestrain = wrap_te(Mesh.PMLElastic2d_te,   'pstrain')
Mesh.PMLElastic2d_te_planestress = wrap_te(Mesh.PMLElastic2d_te,   'pstress')
Mesh.PMLElasticAxis_te           = wrap_te(Mesh.PMLElasticAxis_te, 'axis')
Mesh.PMLELastic3d_te             = wrap_te(Mesh.PMLElastic3d_te,   '3d')


-- Set ordinary elastic materials
-- 
local function wrap_m(constructor, etype)
  local compute_D = {
    ['pstrain'] = elasticity_2D_strain,
    ['pstress'] = elasticity_2D_stress,
    ['axis']    = elasticity_axis,
    ['3D']      = elasticity_3D
  }

  local compute_D_cubic = {
    ['pstrain'] = cubic_elasticity_2D_strain,
    ['pstress'] = cubic_elasticity_2D_stress,
    ['3d']      = cubic_elasticity_3D
  }

  return function(self, mtype, wafer, angle)  
    mtype = fill_mech(get_material(mtype))
    local D = {}
    if wafer and angle then
      local f = compute_D_cubic[etype] or
                error('Cannot make te cubic ' .. etype .. ' material')
      f(D, mtype.lambda, mtype.mu, mtype.alpha, wafer_orientation(wafer,angle))
    else
      local lambda, mu = get_lame(mtype)
      compute_D[etype](D, lambda, mu)
    end
    return constructor(self, D, mtype.rho)
  end
end

Mesh.PMLElastic2d_planestrain = wrap_m(Mesh.PMLElastic2d,   'pstrain')
Mesh.PMLElastic2d_planestress = wrap_m(Mesh.PMLElastic2d,   'pstress')
Mesh.PMLElasticAxis           = wrap_m(Mesh.PMLElasticAxis, 'axis')
Mesh.PMLElastic3d             = wrap_m(Mesh.PMLElastic3d,   '3d')


---------------------------

-- Add an electrode tying the voltage at nodenum to a new global given
-- by efunc.  lt defines a scaling factor for going between the two.
-- 
function add_electrode(mesh,etype,nodenum,efunc,lt)
  local constructors = {
    ['electrode'  ] = mesh.Electrode,
    ['electrode2' ] = mesh.Electrode2
  }
  etype = string.lower(etype)
  assert(constructors[etype], 'Unknown element type ' .. etype)
  local idg1   = mesh:add_global(efunc,'V','Q2')
  local elt    = constructors[etype](mesh,idg1,lt or 1) 
  local eltnum = mesh:add_element({nodenum}, elt)
  return eltnum,idg1  
end


---------------------------

local function add_assembler_bc_wrapper(mesh)
  if mesh.u_bc_list or mesh.f_bc_list then return end
  mesh.u_bc_list = {n = 0}
  mesh.f_bc_list = {n = 0}
  local old_form_bcs = mesh.form_bcs
  function mesh:form_bcs()
    if self.u_bc_list then
      local ubc = QBCAssembler:new(self,'u')
      for i,bcfunc in ipairs(self.u_bc_list) do
        bcfunc(self,ubc)
      end
      ubc:delete()
    end
    if self.f_bc_list then
      local fbc = QBCAssembler:new(self,'f')
      for i,bcfunc in ipairs(self.f_bc_list) do
        bcfunc(self,fbc)
      end
      fbc:delete()
    end
    old_form_bcs(mesh)
  end
end

function add_u_bc(bcfunc)
  add_assembler_bc_wrapper(mesh)
  mesh.u_bc_list.n = mesh.u_bc_list.n + 1
  mesh.u_bc_list[mesh.u_bc_list.n] = bcfunc
end

function add_f_bc(bcfunc)
  add_assembler_bc_wrapper(mesh)
  mesh.f_bc_list.n = mesh.f_bc_list.n + 1
  mesh.f_bc_list[mesh.f_bc_list.n] = bcfunc
end

function clamp_boundary(func, ...)
  add_u_bc(function(mesh,ubc)
    for j = 1,mesh:numnp() do
      if func(mesh:x(j-1)) then
        for k,slot in ipairs(arg) do
          ubc:set(mesh:inode(var_slots[slot],j-1), 0)
        end
      end
    end
  end)
end

function point_load(x,y,p)
  add_f_bc(function(mesh,fbc)
    for j = 1,mesh:numnp() do
      if mesheq(mesh:x(0,j-1),x) and mesheq(mesh:x(1,j-1),y) then
        for key,val in pairs(p) do
          if var_slots[key] then
            fbc:set(mesh:inode(var_slots[key],j-1),val)
          end
        end
      end
    end
  end)
end


---------------------------

function global_indicator(id,val)
  assert(id, 'Must define global number')
  val = val or 1
  return function(mesh,v)
    v:set(mesh:iglobal(id), val)
  end
end

function branch_indicator(dof,id,val)
  assert(dof and id, 'Must define dof and branch number')
  val = val or 1
  return function(mesh,v)
    v:set(mesh:ibranch(dof,id), val)
  end
end

function nodal_indicator(dof,id,val)
  assert(dof and id, 'Must define dof and node number')
  dof = var_slots[dof] or dof
  val = val or 1
  return function(mesh,v)
    v:set(mesh:inode(dof,id), val)
  end
end

function nodal2d_indicator(dof,x,y,val)
  assert(dof and x and y, 'Must define dof and node coordinates')
  dof = var_slots[dof] or dof
  val = val or 1
  return function(mesh,v)
    for j = 1,mesh:numnp() do
      if mesheq(x,mesh:x(0,j-1)) and mesheq(y,mesh:x(1,j-1)) then
        v:set(mesh:inode(dof,j-1), val)
      end
    end
  end
end
