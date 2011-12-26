% Driver to compare simulated frequencies and Q values with
% measured values reported in the Michigan Transducers 03 paper
% on diamond disk.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_sweep0a.m,v 1.1 2006/05/01 05:06:55 tkoyama Exp $

param = []; 
param.order = 3;
param.dense = 1e-6/4;
w0 = 2*pi*47.5e6;
%t = linspace(1.2e-6, 1.8e-6, 61);
t1 = linspace(1.2e-6,    1.5e-6, 31);
t2 = linspace(1.505e-6,  1.6e-6, 20);
t3 = linspace(1.61e-6,   1.8e-6, 20);
t  = [t1, t2, t3];

[w,Q] = driver_sweep(w0, 2, 'hdisk', t, param);

w = w / 2e6 / pi;  % Convert to MHz

% -- Sort by Q -- 
for k = 1:length(t)
  [Q2k,I] = sort(Q(:,k));
  Q(:,k) = Q2k;
  w(:,k) = w(I,k);
end

% -- Plot Q --   
figure(1); clf;
semilogy(t/1e-6, [max(Q); min(Q)]);
xlabel('Film thickness (um)');
ylabel('Q');

% -- Plot eigenvalues --
figure(2);  clf;
plot( real(w(1,:)), imag(w(1,:)), 'b.', ...
      real(w(2,:)), imag(w(2,:)), 'g.');
xlabel('Imaginary part (MHz)');
ylabel('Real part (MHz)');

