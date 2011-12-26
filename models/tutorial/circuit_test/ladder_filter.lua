-- Test File for a ladder filter

require 'common.lua'
require 'circuit.lua'

-- Define physical nodes
x1 = {0.0, 0.0,   -- Node 0
      0.0, 2.0,   -- Node 1
      0.0, 3.0,   -- Node 2
      1.0, 2.0,   -- Node 3
      1.0, 3.0,   -- Node 4
      2.0, 3.0,   -- Node 5
      3.0, 3.0,   -- Node 6
      4.0, 2.0,   -- Node 7
      4.0, 3.0,   -- Node 8
      5.0, 0.0,   -- Node 9
      5.0, 3.0,   -- Node 10
      6.0, 0.0,   -- Node 11
      6.0, 1.0,   -- Node 12
      6.0, 2.0,   -- Node 13
      6.0, 3.0,   -- Node 14
      7.0, 0.0,   -- Node 15
      7.0, 3.0}   -- Node 16

-- Construct mesh and add nodes
mesh   = Mesh:new(2)
numnp  = table.getn(x1)/2
mesh:add_node(x1,numnp)

-- Define circuit elements
insert1 = true
insert2 = true

-- LCR-C(Resonator 1)
L1  = 1 
R1  = 0.01
C1  = 1
C1a = 30

-- LCR-C(Resonator 2)
L2  = 1.03
R2  = 0.01
C2  = 1
C2a = 30

-- Source, Load Resistors and Voltage
Rs = 0.001
Rl = 100.0
E  = 1

-- Add circuit elements to mesh
eno1 = circuit.wire{ 0,1 }
circuit.R{ 1,2; R = Rs }

-- Insert resonator 1
if insert1 then
  circuit.wire{ 3,4 }
  circuit.wire{ 7,8 }
  circuit.C{ 3,7; C=C1a  }
  circuit.L{ 4,5; L=L1   }
  circuit.R{ 5,6; R=R1   }
  circuit.C{ 6,8; C=C1   }
else
  circuit.wire{ 4,5 }
  circuit.wire{ 5,6 }
  circuit.wire{ 6,8 }
end

-- Insert resonator 2
if insert2 then
  circuit.C{  9,10; C=C2a }
  circuit.C{ 11,12; C=C2  }
  circuit.R{ 12,13; R=R2  }
  circuit.L{ 13,14; L=L2  }
end

circuit.R{ 15,16; R=Rl }

circuit.wire{  0, 9 }
circuit.wire{  9,11 }
circuit.wire{ 11,15 }

circuit.wire{  2, 4 }
circuit.wire{  8,10 }

circuit.wire{ 10,14 }
circuit.wire{ 14,16 }

circuit.ground{ 0 }

-- Define drive and sense functions
ac_source   = circuit.current_probe(eno1)
bode_senseA = circuit.voltage_probe(1)
bode_senseB = circuit.voltage_probe(16)

