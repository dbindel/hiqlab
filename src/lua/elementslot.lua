-- Var slot management

var_slots = {nslots = 0}

function get_var_slot(name, s1, s2)

  -- Argument checks
  if var_slots[name] then return var_slots[name] end
  assert(s1, 'Variable scale must be specified')
  assert(s2, 'Dual variable scale must be specified')

  -- Ensure that the dim_scales.vars array is set up
  dim_scales      = dim_scales      or {}
  dim_scales.vars = dim_scales.vars or {}

  -- Allocate a new slot and add to vars
  local slot = var_slots.nslots
  var_slots.nslots = var_slots.nslots + 1
  var_slots[name] = slot
  dim_scales.vars[2*slot+1] = s1
  dim_scales.vars[2*slot+2] = s2

  return slot

end

function Element:set_slots(...)
  local j = 0
  for i = 1,arg.n,3 do
    self:id_slot(j, get_var_slot(arg[i], arg[i+1], arg[i+2]))
    j = j + 1
  end
end


-- Wrappers around standard constructors

local function wrapslots(p)
  local etype  = p[1]
  local nodal  = p.nodal
  local branch = p.branch
  local constructor = Mesh[etype]
  Mesh[etype] = function(...)
    local elt = constructor(unpack(arg))
    local j = 0
    for i,v in ipairs(nodal) do
      for k = 1,table.getn(v),3 do
        elt:id_slot(j, get_var_slot(v[k], v[k+1], v[k+2]))
        j = j + 1
      end
    end
    elt.aux_scales = branch
    return elt
  end
end

-- Electrostatics = scalar with zero mass term
Mesh.Electrostatic1d   = function(m,k) return m:PMLScalar1d(k,0)   end
Mesh.Electrostatic2d   = function(m,k) return m:PMLScalar2d(k,0)   end
Mesh.ElectrostaticAxis = function(m,k) return m:PMLScalarAxis(k,0) end
Mesh.Electrostatic3d   = function(m,k) return m:PMLScalar3d(k,0)   end

local mech2d_slots  = {'ux', 'L', 'F2', 'uy', 'L', 'F2'}
local mech2hd_slots = {'ux', 'L', 'F',  'uy', 'L', 'F'}
local mech3d_slots  = {'ux', 'L', 'F',  'uy', 'L', 'F',  'uz', 'L', 'F'}
local circuit_slots = {'vol', 'V', 'A'}
local phi1d_slots   = {'phi', 'V', 'Q1'}
local phi2d_slots   = {'phi', 'V', 'Q2'}
local phi3d_slots   = {'phi', 'V', 'Q'}
local temp2d_slots  = {'theta', 'Th', 'Qt2'}
local temp3d_slots  = {'theta', 'Th', 'Qt'}

wrapslots{'Electrode',  nodal = {circuit_slots}, branch = {'Q2','V'}}
wrapslots{'Electrode2', nodal = {circuit_slots}, branch = {'Q2','V','A','V'}}
wrapslots{'Resistor',   nodal = {circuit_slots}}
wrapslots{'Capacitor',  nodal = {circuit_slots}}
wrapslots{'Capacitor1', nodal = {circuit_slots}}
wrapslots{'Capacitor2', nodal = {circuit_slots}}
wrapslots{'Inductor',   nodal = {circuit_slots}, branch = {'A','V'}}
wrapslots{'VIsrc',      nodal = {circuit_slots}, branch = {'A','V'}}

wrapslots{'CoupleEM2d', nodal = {mech2d_slots, phi2d_slots}}

wrapslots{'Electrostatic1d',   nodal = {phi1d_slots}}
wrapslots{'Electrostatic2d',   nodal = {phi2d_slots}}
wrapslots{'ElectrostaticAxis', nodal = {phi3d_slots}}
wrapslots{'Electrostatic3d',   nodal = {phi3d_slots}}

wrapslots{'Elastic2d',         nodal = {mech2d_slots}}
wrapslots{'PMLElastic2d',      nodal = {mech2d_slots}}
wrapslots{'PMLElasticAxis',    nodal = {mech2hd_slots}}
wrapslots{'PMLElasticTAxis',   nodal = {mech2hd_slots}}
wrapslots{'PMLElastic3d',      nodal = {mech3d_slots}}

wrapslots{'PMLElastic2d_pz',   nodal = {mech2d_slots,  phi2d_slots}}
wrapslots{'PMLElastic2hd_pz',  nodal = {mech2hd_slots, phi3d_slots}}
wrapslots{'PMLElastic3d_pz',   nodal = {mech3d_slots,  phi3d_slots}}

wrapslots{'PMLElastic2d_te',   nodal = {mech2d_slots,  temp2d_slots}}
wrapslots{'PMLElasticAxis_te', nodal = {mech2hd_slots, temp3d_slots}}
wrapslots{'PMLElastic3d_te',   nodal = {mech3d_slots,  temp3d_slots}}

