% Driver to compare simulated frequencies and Q values with
% measured values reported in the Michigan Transducers 03 paper
% on diamond disk.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_sweep0.m,v 1.1 2006/05/01 05:04:09 tkoyama Exp $
 
param.disk_material = 'diamond';
param.rpost =  1.6e-6 / 2;
param.rdisk = 24.0e-6 / 2;
param.hdisk =  3e-6;
param.hpost =  8e-7;
param.order = 3;
param.dense = 1e-6/6;
w0 = 2*pi*455e6;
t = linspace(2.4e-6, 2.6e-6, 101);

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

% -- Zoomed plot of eigenvalues --
figure(3);  clf; hold on;
I = [31, 36, 41, 46, 51];
hx = xlabel('Real frequency (MHz)');
hy = ylabel('Imag frequency (MHz)');
set(hx, 'FontSize', 12, 'FontWeight', 'bold');
set(hy, 'FontSize', 12, 'FontWeight', 'bold');
plot(w.');
plot(w(:,I).', '*');
axis([452.5 454.5 0 1.6]);
ls = {};
for k = 1:length(I)
  ls{k} = sprintf('(%c) %.2f um', 'a'-1+k, t(I(k))*1e6);
  coord = w(:,I(k));
  xx = real(coord)-0.05;
  yy = imag(coord);
  h = text(xx, yy, sprintf('%c', 'a'-1+k));
  set(h, 'FontSize', 12, 'FontWeight', 'bold');
end
h = text(454, 0.4, ls);
set(h, 'FontSize', 12, 'FontWeight', 'bold');
