function [Ha] = ladder_filter_analytical(L,freq);

% -- Compute analytical result
Rs  = Lua_get_double(L,'Rs');
Rl  = Lua_get_double(L,'Rl');
pL1 = Lua_get_double(L,'L1');
pC1 = Lua_get_double(L,'C1');
pC1a= Lua_get_double(L,'C1a');
pR1 = Lua_get_double(L,'R1');
pL2 = Lua_get_double(L,'L2');
pC2 = Lua_get_double(L,'C2');
pC2a= Lua_get_double(L,'C2a');
pR2 = Lua_get_double(L,'R2');

Ha  = zeros(1,size(freq,2));
for i=1:size(freq,2)
    w  = freq(i);

    Z1a= 1/complex(0,w)/pC1 + pR1 + complex(0,w)*pL1;
    Z1b= 1/complex(0,w)/pC1a;
    Z2a= 1/complex(0,w)/pC2 + pR2 + complex(0,w)*pL2;
    Z2b= 1/complex(0,w)/pC2a;

    Z1 = Z1a*Z1b/(Z1a+Z1b);
    Z2 = Z2a*Z2b/(Z2a+Z2b);
    Z3 = Z2*Rl/(Z2+Rl);
    Zt = Rs + Z1 + Z3;

    Ha(i) = (Rs+Rl)/Rl*Z3/Zt;
end;
