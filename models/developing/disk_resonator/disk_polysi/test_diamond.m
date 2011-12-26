 
% Driver to compare simulated frequencies and Q values with
% measured values reported in the Michigan Transducers 03 paper
% on diamond disk.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_diamond.m,v 1.1 2006/05/01 05:04:09 tkoyama Exp $
 
fid = fopen('test_diamond.txt', 'w');

DIAMOND = 0;
SILICON = 1;

device_info = ...
  [DIAMOND, 1.6, 24,  455.8,  24117;   % 1
   DIAMOND, 1.6, 24, 1272.2,  12050;   % 2
   SILICON, 1.6, 22,  245.1,   8100;   % 3
   DIAMOND, 1.6, 22,  497.58, 55300;   % 4
   SILICON, 1.6, 22,  657.45,  6400;   % 5
   DIAMOND, 1.6, 22, 1388.1,  10680;   % 6
   SILICON, 2.0, 22,  243.8,   4500;   % 7
   DIAMOND, 2.0, 22,  495.64, 27500;   % 8
   SILICON, 2.0, 22,  655.82,  4900;   % 9
   DIAMOND, 2.0, 22, 1386.9,   8760;   % 10
   DIAMOND, 1.6, 20,  545.9,  17458;   % 11
   DIAMOND, 1.6, 20, 1507.8,  11555;   % 12
   DIAMOND, 2.0, 20,  547.2,  11448;   % 13
   DIAMOND, 2.0, 20, 1519.7,   4648;   % 14
   DIAMOND, 1.7, 17,  639.2,  10531;   % 15
   DIAMOND, 2.0, 16,  681.7,  10246];  % 16

I = 1:size(device_info,1);
[tmp, J] = sort(device_info(I,4));  I = I(J);
[tmp, J] = sort(device_info(I,1));  I = I(J);

hlist = linspace(2e-6, 3e-6, 11);

param.hpost = 8e-7;

fprintf(fid, 'Material  stem  disk  thick f_meas    f_comp  Q_meas  Q_comp\n'); 
fprintf(fid, '--------  ----  ----  ----- ------    ------  ------  ------\n');

ntotal = length(device_info)*length(hlist);
ndone = 0;

t0 = clock;
for kk = 1:length(device_info)
  
  k = I(kk);
  if device_info(k,1) == DIAMOND
    param.disk_material = 'diamond';
  else
    param.disk_material = 'silicon';
  end
  param.order = 3;
  param.dense = 1e-6/6;
  param.rpost = device_info(k,2)*5e-7;
  param.rdisk = device_info(k,3)*5e-7;
  w0          = device_info(k,4)*2e6*pi;

  for ll = 1:length(hlist)
 
    param.hdisk = hlist(ll); 
    [w,Q] = driver_mode(w0, param);
    device_comp(k,1) = real(w/2e6/pi);
    device_comp(k,2) = imag(w/2e6/pi);
    device_comp(k,3) = Q;
    ndone = ndone + 1;
    t1 = clock;
    te = etime(t1,t0);
  
    fprintf('Material: %s\n', param.disk_material);
    fprintf('Diameter: %g %g\n', device_info(k,2), device_info(k,3));
    fprintf('Thickness: %g\n', param.hdisk);
    fprintf('Frequency: %g\t/\t%g\n', w0/2e6/pi, real(w)/2e6/pi);
    fprintf('Q        : %g\t/\t%g\n', device_info(k,5), Q);
    fprintf('Progress : %g%%\n', ndone / ntotal * 100);
    fprintf('Time     : %g elapsed + %g remaining = %g total\n', ...
           te, te*(ntotal-ndone)/ndone, te*ntotal/ndone);

    fprintf(fid, ' %s  ', param.disk_material);
    fprintf(fid, '% 3.1f  ', device_info(k,2));
    fprintf(fid, '% 4.1f  ', device_info(k,3));
    fprintf(fid, '% 3.1f  ', param.hdisk * 1e6);
    fprintf(fid, '% 8.2f  ', device_info(k,4));
    fprintf(fid, '% 8.2f  ', device_comp(k,1));
    fprintf(fid, '% 6.0f  ', device_info(k,5));
    fprintf(fid, '% 6.0f  ', device_comp(k,3));
    fprintf(fid, '\n');

  end

end

save -ascii test_diamond.mat device_info device_comp
fclose(fid);
