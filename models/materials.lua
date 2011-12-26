-- HiQLab
-- Copyright (c): Regents of the University of California

function fill_mech(p)

  -- Fill in crystal property
  p.crystal = p.crystal or 'isotropic'

  if (p.crystal == 'isotropic') or (p.crystal == 'cubic') then
    -- Fill in E and nu from mu, lambda
    p.E  = p.E  or p.mu*(3*p.lambda+2*p.mu)/(p.lambda+p.mu)
    p.nu = p.nu or p.lambda/2.0/(p.lambda+p.mu)

    -- Fill in lambda, mu, and alpha from E and nu
    p.lambda =  p.lambda or p.E*p.nu/(1.0 - 2.0*p.nu)/(1.0 + p.nu)
    p.mu     =  p.mu     or p.E/2.0/(1.0 + p.nu)
    p.alpha  =  p.alpha  or p.lambda + 2*p.mu

    -- Fill equivalent hexagonal properties
    p.c11    = p.c11 or (p.lambda + 2*p.mu)
    p.c12    = p.c12 or  p.lambda
    p.c13    = p.c13 or  p.lambda
    p.c33    = p.c33 or (p.lambda + 2*p.mu)
    p.c55    = p.c55 or  p.mu

  elseif p.crystal == 'hexagonal' then
    p.E      = p.E or p.c33
  else
    error('Unknown crystal type: ' .. p.crystal)
  end

  p.rho = p.rho or 0

  return p
end


function fill_piezo(p)

  -- Fill in crystal property
  p.crystal = p.crystal or 'isotropic'
  p.kds1    = p.kds1    or 0
  p.kds3    = p.kds3    or 0
  p.d16     = p.d16     or 0
  p.d31     = p.d31     or 0
  p.d33     = p.d33     or 0

  if (p.crystal == 'isotropic') or (p.crystal == 'cubic') then
    p.d   = p.d33  or 0     -- Maybe should change this
    p.kds = p.kds3 or 0
  elseif p.crystal == 'hexagonal' then
    -- Fill e16,e31,e33 from d16,d31,d33
    p.e16 =  p.e16 or   p.c55 * p.d16
    p.e31 =  p.e31 or  (p.c11 + p.c12) * p.d31 + p.c13 * p.d33
    p.e33 =  p.e33 or 2*p.c13 * p.d31          + p.c33 * p.d33

    -- Fill d16,d31,d33 from e16,e31,e33
    local det = (p.c11+p.c12)*p.c33 - 2*p.c13*p.c13
    p.d16 = p.d16 or     p.e16/p.c55
    p.d31 = p.d31 or (   p.e31*p.c33 - p.e33*p.c13)/det
    p.d33 = p.d33 or (-2*p.e31*p.c13 + p.e33*(p.c11+p.c12))/det

    p.d   = p.d33
    p.kds = p.kds3
  else
    error('Unknown crystal type: ' .. p.crystal)
  end

  return p
end


-- Set table for non-dimensionalization constants
dim_scales      = dim_scales or {}

-- Set ambient temperature and permittivity of free space
dim_scales.T0   = 293.15
dim_scales.eps  = 8.85419e-12


-------------------------------------------------
--           Material Database
-------------------------------------------------

sc_silicon = fill_mech({
  -- SC Silicon properties
  crystal= 'cubic',
  rho    = 2330.0,
  lambda =   64.0e9,
  mu     =   80.0e9,
  alpha  =  166.0e9, 
  at     =    2.6e-6,
  cp     =  712.0,
  kt     =  148.0
})


poly_silicon = fill_mech({
  E      =  150.0e9,
  nu     =  0.226,
  rho    = 2330.0,
  at     =   2.6e-6,
  cp     =  712.0,
  kt     =   30.0
})


amor_silicon = fill_mech({
  -- Pseudo amorphous Silicon parameters
  -- Have taken sc_silicon and modified kt
  crystal= 'cubic',
  rho    = 2330.0,
  lambda =   64.0e9,
  mu     =   80.0e9,
  alpha  =  166.0e9,
  at     =   2.6e-6,
  cp     =  712.0,
  kt     =    1.0,
})


silicon = fill_mech({
  -- Polysilicon properties taken from Michigan Si disk paper
  rho = 2330,
  E   = 150e9,
  nu  = 0.3
})


silicon2 = fill_mech({
  -- Polysilicon properties (local)
 rho    = 2300,
 E      = 165e9,
 nu     = 0.3,
 at     = 2.6e-6,
 cp     =  712.0,
 kt     =   30.0
})


hfo2 = fill_mech({
  -- Hafnium Oxide properties(hypothetical, should correct)
  rho = 9680,
  E   = 240e9,
  nu  = 0.3
})


sic = fill_mech({
  -- Polysilicon carbide properties (from Sunil -- SCS values?)
 rho = 3100,
 E   = 700e9,
 nu  = 0.25,
})


sige = fill_mech({
  -- Estimated properties for Si_0.4 Ge_0.6 film
  -- c = 5800 m/s (about 6000?)
  -- E   = 139e9; -- Original estimate
  rho = 4125;
  E   = 139e9;
  nu  = 0.28;
})


siox = fill_mech({
  -- SiO2 from http://www.memsnet.org/material/silicondioxidesio2film/
  rho = 2200,
  E   = 70e9,
  nu  = 0.17
})


diamond = fill_mech({
  -- Polydiamond paper from Journal of Applied Physics
  rho = 3470,
  E   = 1120e9,
  nu  = 0.09
})


aln = {
  -- Aluminum Nitride Properties from ?????
  -- Actually a hexagonal structure 5 elastic constants
  -- (??)
  -- c11,c12,c13,c33,c44, (c11-c12)/2
  -- [c11,c12,c13,          0,  0,  0,
  --     ,c11,c13,          0,  0,  0,
  --         ,c33,          0,  0,  0,
  --             ,(c11-c12)/2,  0,  0,
  --                         ,c55,  0,
  --                             ,c55]
  crystal = 'hexagonal',
  rho = 3260,
  c11 = 345e9,
  c12 = 125e9,
  c13 = 120e9,
  c33 = 395e9,
  c55 = 118e9,
--  c11 = 388.957816377171,
--  c12 = 122.828784119107,
--  c13 = c12,
--  c33 = c11,
--  c55 = 133.064516129032,
  -- 3 piezoelectric constants
  -- Differs in bulk and poly state(current values are bulk)
  -- (??) Still don't which values to use
  -- [   0,   0,   0,  0,   0,  d16,
  --     0,   0,   0,  0, d16,    0,
  --   d31, d31, d33,  0,   0,    0]
  piezo  = 1,
  d16    =-4.07e-12,
  d33    = 5.53e-12,
  d31    =-2.65e-12,
  -- 2 dielectric constants at constant stress(static)
  --                           at high freq 4.6
  -- [ kds1,    0,     0,
  --      0, kds1,     0,
  --      0,    0,  kds3]
  kds1   = 8.5 * 8.85419e-12,
  kds3   = 4.5 * 8.85419e-12
}
aln = fill_mech(aln)
aln = fill_piezo(aln)


aln_isotropic = {
  -- Mechanically isotropic but has piezoelectricity
  -- pseudo material
  crystal= 'isotropic',
  rho    = 3260,
  E      = 330,
  nu     = 0.24,
  piezo  = 1,
  d16    =-4.07e-12,
  d31    =-2.65e-12,
  d33    = 5.53e-12,
  kds1   = 8.5 * 8.85419e-12,
  kds3   = 8.5 * 8.85419e-12
}
aln_isotropic = fill_mech(aln_isotropic)
aln_isotropic = fill_piezo(aln_isotropic)


aln_piazza = {
  crystal= 'hexagonal',
  rho    = 3200,
  c11    = 410e9,
  c12    = 140e9,
  c13    = 100e9,
  c33    = 390e9,
  c55    = 120e9,
  e16    =-0.5,
  e31    =-0.75,
  e33    = 1,
  kds1   = 9 * 8.85418e-12,
  kds3   = 9 * 8.85418e-12
}
aln_piazza = fill_mech (aln_piazza)
aln_piazza = fill_piezo(aln_piazza)


pt_piazza = {
  crystal= 'isotropic',
  rho    = 21090,
  E      = 168e9,
  nu     = 0.38,
  resist = 2.4e-7  -- Ohm m
}
pt_piazza = fill_mech (pt_piazza)
pt_piazza = fill_piezo(pt_piazza)


al_piazza = {
  crystal= 'isotropic',
  rho    = 2700,
  E      = 70e9,
  nu     = 0.35,
  resist = 5e-8   -- Ohm m
}
al_piazza = fill_mech (al_piazza)
al_piazza = fill_piezo(al_piazza)
