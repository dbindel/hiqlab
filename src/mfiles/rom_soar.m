% [Mk,Dk,Kk,Lk,Bk,Vk] = rom_soar(M,D,K,L,B,kk,w0,opt)
%
% Build Second Order Krylov Subspace from matrices
%
% M0, D0, K0
%
% Example:TED case
%
% M0 =  [ Muu    0 ]  D0 = [ 0      0 ]  K0 = [ Kuu Kut]
%       [  0     0 ]       [ Ctu  Ctt ]       [  0  Ktt]
%
% Kh = -w0^2 * M0 +     iw0 * D0 + K0
% Dh =            - 2 * iw0 * M0 + D0
% Mh =                             M0
% A  = -Kh^{-1} * Dh
% B  = -Kh^{-1} * Mh
% Do kk steps of SOAR for the shift-and-inverted eigenvalue
% problem starting from B and project everything onto the
% resulting space.
%
% r_0 = b
% r_1 = A * r0
%     = - Kh^{-1} * Dh * r0
% r_j = A * r_(j-1) + B * r_(j-2)
%     = - Kh^{-1} * ( Dh * r_(j-1) + Mh * r_(j-2)  )
%
% Options:
%   realbasis (0)    - Transform to a real-valued basis
%   realtol   (1e-8) - When converting to real basis, discard
%                      directions with singular values < realtol
%   structurep(0)    - structure preserving if 0
%   mechdof(0)       - Free dofs for mechanical portion
%                      this field is required for structurep opt
function [Mk,Dk,Kk,Lk,Bk,Vk] = rom_soar(M,D,K,L,B, kk, w0, opt)

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 6, w0 = 0;   end
if nargin < 7, opt = []; end

realbasis = qoptdefault(opt, 'realbasis', 0   );
realtol   = qoptdefault(opt, 'realtol',   1e-8);
structurep= qoptdefault(opt, 'structurep',0   );
nm        = qoptdefault(opt, 'mechdof',   0   );
if nm==0
   structurep = 0;
end;

if kk <= 0
  Mk = M;
  Dk = D;
  Kk = K;
  Lk = L;
  Bk = B;
  Vk = [];
  return;
end

[L0,U0,P0,Q0] = qlu(K +   complex(0,w0) * D - w0^2*M);

n = length(M);
H = zeros(kk+1,kk);
V = zeros(n,kk+1);

if ~isempty(B)
  r= Q0*(U0\(L0\(P0*B)));
else
  r= rand(n,1);
end
r= r/norm(r);

V(:,1) = r;
f      = zeros(n,1);

for k = 1:kk-1
  r = -Q0*(U0\(L0\(P0*(  M*f + D*V(:,k) + 2*complex(0,w0)*M*V(:,k)  ))));
  for i = 1:k
    H(i,k) = V(:,i)'*r;
    r = r-H(i,k)*V(:,i);
  end
  H(k+1,k) = norm(r);
  V(:,k+1) = r/H(k+1,k);
  e        = zeros(k,1);
  e(k,1)   = 1;
  f        = V(:,1:k)*(H(2:k+1,1:k)\e);
end

% -- Define projection basis
if     ~realbasis & ~structurep
  Vk = V(:,1:kk);
elseif ~realbasis &  structurep
  VVm= V(1:nm  ,1:kk);
  [U1,S1,V1] = svd(VVm,0);
  S1 = diag(S1);
  Vkm= U1(:,find(S1 > realtol));

  VVt= V(nm+1:n,1:kk);
  [U1,S1,V1] = svd(VVt,0);
  S1 = diag(S1);
  Vkt= U1(:,find(S1 > realtol));

  Vk = [              Vkm,             zeros(size(Vkm,1),size(Vkt,2));
        zeros(size(Vkt,1),size(Vkm,2)),                       Vkt    ;];
elseif  realbasis & ~structurep
  VV = [real(V(:,1:kk)), imag(V(:,1:kk))];
  [U1,S1,V1] = svd(VV,0);
  S1 = diag(S1);
  Vk = U1(:,find(S1 > realtol));
elseif realbasis  &  structurep
  VV = [real(V(:,1:kk)), imag(V(:,1:kk))];

  VVm= VV(1:nm  ,:);
  [U1,S1,V1] = svd(VVm,0);
  S1 = diag(S1);
  Vkm= U1(:,find(S1 > realtol));

  VVt= VV(nm+1:n,:);
  [U1,S1,V1] = svd(VVt,0);
  S1 = diag(S1);
  Vkt= U1(:,find(S1 > realtol));

  Vk = [              Vkm,             zeros(size(Vkm,1),size(Vkt,2));
        zeros(size(Vkt,1),size(Vkm,2)),                       Vkt    ;];
end

% -- Project onto the Arnoldi basis
Mk = Vk'*M*Vk;
Dk = Vk'*D*Vk;
Kk = Vk'*K*Vk;
if isempty(L), Lk = []; else Lk = Vk'*L; end
if isempty(B), Bk = []; else Bk = Vk'*B; end

