% [L,U,P,Q] = qlu(A,opt)
%
% Calls UMFPACK if it exists. If not the
% R12 LU or the R13-R14 LU to factor a sparse A
% with both row and column pivoting.
% Input:
%   A      - Sparse matrix
%
% Output:
%   L,U,P,Q - sparse matrices LU = PAQ
%
% Options are
%   use_umfpack - Use UMFPACK if present? (default 1)
%
% Other options are passed directly to eigs (if it's used).  If
% disp is not explicitly set, it is assumed to be zero (not the
% more verbose default of 1 generally used by eigs).

function [L,U,P,Q] = qlu(A,opt)

% HiQLab
% Copyright (c): Regents of the University of California
if nargin < 2, opt = []; end;
use_umfpack = qoptdefault(opt, 'use_umfpack', 1);

if exist('umfpack') == 3 & use_umfpack
  [L,U,P,Q] = umfpack(A);
elseif exist('splu')
  [L,U,P,Q] = splu(A);
elseif str2num(version('-release')) >= 13
  [L,U,P,Q] = lu(A);
else
  n  = size(A,1);
  iq = colamd(A);
  i1 = ones(n,1);
  [L,U,P] = lu(A(:,iq));
  Q = sparse(iq,1:n,i1, n,n);
end

