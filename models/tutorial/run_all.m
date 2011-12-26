cd disk
circular_disk
cd ..

cd patch_test
patch_test
cd ..

cd cantilever
cant_te_sta
cant_te_dyn
cant_wa_m_tra
cant_wa_m_sta
cant_wa_m_dyn
cant_m_tra
cant_m_sta
cant_wa_te_tra
cant_m_dyn
cant_wa_te_sta
cant_em_tra
cant_em_sta
cant_wa_te_dyn
%cant_em_dyn
cant_te_tra
cd ..

cd circuit_test
low_pass_filter
ladder_filter
cd ..

% ./beam_sample/beamsweep.m
% ./beam_sample/beammodes.m

cd pml1d
%pml1d
cd ..

cd arch
arch
cd ..

cd pmlaxis
pmlaxis
pmlaxis_s
cd ..

cd pml2d
pml2d
pml2ds
cd ..

% ./circular_disk/circular_disk.m
%./plate2d/beamsweep.m
%./plate2d/beammodes.m

qcloseall;
