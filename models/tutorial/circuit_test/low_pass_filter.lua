require 'common.lua'
require 'circuit.lua'

R = 3;  -- Resistance
C = 7;  -- Capacitance
E = 1;  -- AC voltage magnitude

-- Construct mesh and add nodes
mesh = Mesh:new(2);
idx1 = mesh:add_node{0,0}
idx2 = mesh:add_node{1,0}
idx3 = mesh:add_node{1,1}
idx4 = mesh:add_node{0,1}

-- Add elements to mesh
circuit.wire{ idx1,idx2 }
circuit.C{ idx2,idx3; C=C }
circuit.R{ idx3,idx4; R=R }
elem1 = circuit.wire{ idx4,idx1 }
circuit.ground{ idx1 }

-- Current source and sense probes
ac_source   = circuit.current_probe(elem1, E)
bode_senseA = circuit.voltage_probe(idx3)
bode_senseB = circuit.voltage_probe(idx4)

