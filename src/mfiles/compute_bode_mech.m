% [H] = compute_bode_mech(freq, M,K,L,B, opt)
%
% Compute and optionally plot a Bode plot of the given system.
% Use the transfer function
%   H(s) = L^T (K-s^2 M)^(-1) B
% evaluated at points in freq.
%
% Options are:
%   k      ([]) - Dimension of reduced model to be used (default to none)
%   w0     (0)  - Shift frequency used in model reduction
%   plot   (0)  - Should the result be plotted?
%   usehz  (0)  - Assume input frequencies are Hz (else in radian/s)
%   schur  (1)  - Transform to Schur form for faster computation (only dense)
%   verb   ([]) - Print status as we compute (default: true if sparse)
%
% Other options can be passed through opt to either plot_bode or rom_arnoldi.
% TODO: How should the transfer function be automatically re-dimensionalized?

function [H] = compute_bode_mech(freq, M,K,L,B, opt)

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 6, opt = []; end

k   = qoptdefault(opt, 'k',      []);
w0  = qoptdefault(opt, 'w',      0 );
plt = qoptdefault(opt, 'plot',   0 );
v   = qoptdefault(opt, 'verb',   []);
sch = qoptdefault(opt, 'schur',  1 );
hz  = qoptdefault(opt, 'usehz',  0 );
n   = length(freq);

if hz
  freq = freq * 2*pi;
  w0   = w0   * 2*pi;
  opt.usehz = 0;
end

if isempty(v)
  v = issparse(M);
end

if ~isempty(k)
  [M,K,L,B] = rom_arnoldi(M,K,L,B, kk, w0, opt);
end

if sch & ~issparse(M)
  [U,T] = schur(M\K);
  L = U'*L;
  B = U'*(M\B);
  K = T;
  M = eye(length(T));
end

H = zeros(n,1);

for idx = 1:n
  w = freq(idx);
  if v
    fprintf('%d: %d Hz\n', idx, freq(idx)/2/pi);
  end
  H(idx) = L' * ((K - w^2*M) \ B);
end

if plt
  plot_bode(freq, H, opt);
end
