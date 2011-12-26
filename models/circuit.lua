circuit = {}

function circuit.R(p)
  local nodeA = p[1] or p.A
  local nodeB = p[2] or p.B
  local R  = p[3] or p.R
  local Ri = p[4] or p.Ri or 0
  assert(nodeA and nodeB, 'Nodal points must be defined')
  assert(R and Ri, 'Resistance must be defined')
  return mesh:add_element({nodeA, nodeB}, mesh:Resistor(R,Ri))
end

function circuit.C(p)
  local nodeA = p[1] or p.A
  local nodeB = p[2] or p.B
  local C  = p[3] or p.C
  local Ci = p[4] or p.Ci or 0
  assert(nodeA and nodeB, 'Nodal points must be defined')
  assert(C and Ci, 'Capacitance must be defined')
  return mesh:add_element({nodeA, nodeB}, mesh:Capacitor(C,Ci))
end

function circuit.L(p)
  local nodeA = p[1] or p.A
  local nodeB = p[2] or p.B
  local L  = p[3] or p.L
  local Li = p[4] or p.Li or 0
  assert(nodeA and nodeB, 'Nodal points must be defined')
  assert(L and Li, 'Inductance must be defined')
  return mesh:add_element({nodeA, nodeB}, mesh:Inductor(L,Li))
end

function circuit.wire(p)
  local nodeA = p[1] or p.A
  local nodeB = p[2] or p.B
  assert(nodeA and nodeB, 'Nodal points must be defined')
  circuit.wire_elt = circuit.wire_elt or mesh:VIsrc()
  return mesh:add_element({nodeA, nodeB}, circuit.wire_elt)
end

local function add_bc_wrapper(mesh)
  if not mesh.ground_bc_list then
    mesh.ground_bc_list = {n = 0}
    local old_form_bcs = mesh.form_bcs
    function mesh:form_bcs()
      local slot = var_slots['vol'] or 0
      for i,nodei in ipairs(mesh.ground_bc_list) do
        local k = self:inode(slot,nodei[1])
        self:set_bcode(k, 'u')
        self:set_bv(k, nodei[2])
      end
      old_form_bcs(mesh)
    end
  end
end

function circuit.ground(p)
  add_bc_wrapper(mesh)
  for i,nodei in ipairs(p) do
    mesh.ground_bc_list.n = mesh.ground_bc_list.n + 1
    mesh.ground_bc_list[mesh.ground_bc_list.n] = {nodei, 0}
  end
end

function circuit.set_Vdc(node,Vdc)
  add_bc_wrapper(mesh)
  mesh.ground_bc_list.n = mesh.ground_bc_list.n + 1
  mesh.ground_bc_list[mesh.ground_bc_list.n] = {node, Vdc}
end

function circuit.voltage_probe(node,value)
  value = value or 1
  return function(mesh,v)
    local slot = var_slots['vol'] or 0
    local k = mesh:inode(slot,node)
    v:set(k,value)
  end
end

function circuit.current_probe(wire,value)
  value = value or 1
  return function(mesh,v)
    local k = mesh:ibranch(0,wire)
    v:set(k,value)
  end
end
