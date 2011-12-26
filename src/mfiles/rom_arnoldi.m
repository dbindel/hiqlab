% [Mk,Kk,Lk,Bk,Vk] = rom_arnoldi(M,K,L,B, kk, w0, opt)
%
% Do kk steps of ordinary Arnoldi on K - w0^2 * M starting from B
% and project everything onto the resulting space.
% Options:
%   use_umfpack(0)   - Use UMFPACKs LU
%   realbasis (0)    - Transform to a real-valued basis
%   skewproj (0)     - Use skew projection on a complex basis
%   realtol   (1e-8) - When converting to real basis, discard
%                      directions with singular values < realtol

function [Mk,Kk,Lk,Bk,Vk] = rom_arnoldi(M,K,L,B, kk, w0, opt)

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 6, w0 = 0;   end
if nargin < 7, opt = []; end

qopt.use_umfpack = qoptdefault(opt, 'use_umfpack',  1);
realbasis        = qoptdefault(opt, 'realbasis',    0);
skewproj         = qoptdefault(opt, 'skewproj',     0);
realtol          = qoptdefault(opt, 'realtol',   1e-8);

if kk <= 0
  Mk = M;
  Kk = K;
  Lk = L;
  Bk = B;
  Vk = [];
end

[L0,U0,P0,Q0] = qlu(K-w0^2*M,qopt);

n  = length(M);
H = zeros(kk+1,kk);
V = zeros(n,kk+1);

if ~isempty(B)
  r= Q0*(U0\(L0\(P0*B)));
else
  r = rand(n,1);
end
r= r/norm(r);

V(:,1) = r;
for k = 1:kk-1
  r = Q0*(U0\(L0\(P0*(M*V(:,k)))));
  for i = 1:k
    H(i,k) = V(:,i)'*r;
    r = r-H(i,k)*V(:,i);
  end
  H(k+1,k) = norm(r);
  V(:,k+1) = r/H(k+1,k);
end

if realbasis 

  % Build a real-valued basis from the Arnoldi basis
  VV = [real(V), imag(V)];
  [U1,S1,V1] = svd(VV,0);
  S1 = diag(S1);
  Vk = U1(:,find(S1 > realtol));

  Mk = Vk'*M*Vk;
  Kk = Vk'*K*Vk;
  if isempty(L), Lk = []; else Lk = Vk'*L; end
  if isempty(B), Bk = []; else Bk = Vk'*B; end

elseif skewproj

  % Skew project onto the usual Arnoldi basis
  Vk = V(:,1:kk);
  Mk = Vk.'*M*Vk;
  Kk = Vk.'*K*Vk;
  if isempty(L), Lk = []; else Lk = conj(Vk.'*L); end
  if isempty(B), Bk = []; else Bk = Vk.'*B; end

else

  % Project onto the usual Arnoldi basis
  Vk = V(:,1:kk);
  Mk = Vk'*M*Vk;
  Kk = Vk'*K*Vk;
  if isempty(L), Lk = []; else Lk = Vk'*L; end
  if isempty(B), Bk = []; else Bk = Vk'*B; end

end
