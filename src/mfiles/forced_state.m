% function forced_state(mesh,F,w,opt)
%
% Solve for forced state given forcing vector
%
% Input:
%   mesh     - Mesh object
%   F        - Forcing vector(in nondimensional form)
%   w        - Forcing frequency
%   opt      - Optional parameter structure
%     mkc(0)        - Include damping matrix
%     kmax(1)       - Max number of iterations
%     du_tol(1e-15) - Tolerance for convergence for increment
%     R_tol(1e-15)  - Tolerance for convergence for residual
%
function forced_state(mesh,F,w,opt)

if nargin < 4, opt = [];  end;

opt.use_umfpack = qoptdefault(opt, 'use_umfpack',  exist('umfpack') == 3);
opt.mkc = qoptdefault(opt, 'mkc', 0);
kmax    = qoptdefault(opt, 'kmax', 1);
du_tol  = qoptdefault(opt, 'du_tol', 1e-15);
R_tol   = qoptdefault(opt, 'R_tol' , 1e-15);

if ~opt.mkc
  [M,K, su,sf] = Mesh_assemble_mk_nd(mesh);
  [L0,U0,P0,Q0] = qlu(K-w^2*M,opt);
  numid         = Mesh_get_numid(mesh);
  C             = spalloc(numid,numid,0);
else
  [M,K,C, su,sf] = Mesh_assemble_mkc_nd(mesh);
  [L0,U0,P0,Q0] = qlu(K + 1i*w*C - w^2*M,opt);
end

% -- Conduct a linear solve and iterative refinement
U    = zeros(Mesh_get_numid(mesh), 1);
for k = 1:kmax
  R  = ( (K + 1i*w*C - w^2*M)*U - F )./sf;
  du = Q0*(U0\(L0\(P0*R)));
  U  = U - du.*su;
  fprintf('  %d: %e %e\n', k, norm(du), norm(R));
  if norm(du) < du_tol & norm(R) < R_tol, break; end;
end

Mesh_set_u(mesh,U);
