% [V,w,Q] = pml_mode(M,K,w0,nmodes,opt)
%
% Compute damped mode shape, complex frequency, and associated Q.
%
% Input:
%   M      - Mass matrix
%   K      - Stiffness matrix
%   w0     - Estimated frequency (radians/s)
%   nmodes - Number of modes to compute (default is 1)
%   opt    - Options structure
%
% Output:
%   V  - Mode vector (normalized to unit Euclidean length)
%   w  - Mode frequency (radians/s)
%   Q  - Quality factor
%
% Options are
%   use_matlab  - Use MATLAB's eigs?                 (default 0)
%   use_umfpack - Use UMFPACK with eigs, if present? (default 1)
%   linearized  - Is this a linearized eigs problem? (default 0)
%   isreal      - Is this a real problem?            (default 0)
% Other options are passed directly to eigs (if it's used).  If
% disp is not explicitly set, it is assumed to be zero (not the
% more verbose default of 1 generally used by eigs).

function [V,w,Q] = pml_mode(M,K,w0,nmodes,opt)

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 4, nmodes = 1; end
if nargin < 5, opt = [];   end

opt.use_matlab  = qoptdefault(opt, 'use_matlab',  0);
opt.use_umfpack = qoptdefault(opt, 'use_umfpack', 1);
opt.linearized  = qoptdefault(opt, 'linearized',  0);
opt.isreal      = qoptdefault(opt, 'isreal',      0);
opt.disp        = qoptdefault(opt, 'disp',        0);

if opt.linearized
  w0 = sqrt(w0);
end

if opt.use_matlab
  [L1,U1,P1,Q1] = qlu(K-w0^2*M,opt);
  [V,D] = eigs(@afun, length(M), nmodes, w0^2, opt, L1, U1, P1, Q1, M);
  d = diag(D);
else
  [d,V] = pml_mode_2(K-w0^2*M, M, nmodes);
  d = w0^2 + 1./d(1:nmodes);
  V = V(:,1:nmodes);
end

V = V / diag(sqrt(diag(V'*V)));
w = d;
if ~(opt.linearized) w = sqrt(w); end;
Q = abs(w)./(2*imag(w));


% ---
% Apply y = (K-sigma*M)\(M*x)
%
function [y] = afun(x, L1, U1, P1, Q1, M)
y = Q1*(U1\(L1\(P1*(M*x))));
