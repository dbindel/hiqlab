require 'circuit.lua'

-- Add nodes
yn1=-20e-6
yn2=-10e-6
yn3= 10e-6
yn4= 20e-6

idx1 = mesh:add_node{0,yn1}
idx2 = mesh:add_node{0,yn2}
idx3 = mesh:add_node{0,yn3}
idx4 = mesh:add_node{0,yn4}

-- Define global shape function for plates
function plateg_d(x,y)
   if meshbetween(x, -Cw*2/3, Cw*2/3) and mesheq( y,Ch/2) then
                                  return 0, 0, 1, 0;end
end
function plateg_s(x,y)
   if meshbetween(x, -Cw*2/3, Cw*2/3) and mesheq( y,-Ch/2) then
                                  return 0, 0, 1, 0;end
end

-- Add discrete circuit elements
eno_vd = circuit.wire{ idx3, idx4 }
eno_vs = circuit.wire{ idx2, idx1 }

-- Add Electrode element
eno_d,idg_d = add_electrode('electrode',idx3,plateg_d,
                             get_dim_scale('Lt')/get_dim_scale('L'))
eno_s,idg_s = add_electrode('electrode',idx2,plateg_s,
                             get_dim_scale('Lt')/get_dim_scale('L'))

-- Nodal boundary conditions
circuit.ground{ idx1 }
circuit.set_Vdc( idx4, Vf_dc )

-- Element sensing and driving functions
ac_voltage_source2      = circuit.current_probe(eno_vd)
sense_motional_current2 = circuit.current_probe(eno_vs)
electrode_top_V2        = global_indicator(idg_d)
electrode_top_Q2        = branch_indicator(0,eno_d)
