% function static_state(mesh,opt)
%
% Solve for static state on the stiffness matrix
% If the field nonlinear is not specified only 1
% iteration will be conducted.
%
% Input: mesh  - Mesh object
%       *opt
%         nonlinear('NR')- Conduct non-linear solve
%                       'NR' :Newton-Raphson
%                       'MNR':Modified-Newton-Raphson
%         kmax(20)      - Max number of iterations
%         du_tol(1e-8) -  Tolerance for convergence
%                         for increment
%         R_tol (1e-8) -  Tolerance for convergence
%                         for residual
%         U             - Initial starting vector for solve
%
function static_state(mesh,opt)

if nargin < 2, opt = [];  end;

nonlinear = qoptdefault(opt, 'nonlinear', 'NR');
kmax      = qoptdefault(opt, 'kmax', 20);
du_tol    = qoptdefault(opt, 'du_tol', 1e-8);
R_tol     = qoptdefault(opt, 'R_tol' , 1e-8);
U         = qoptdefault(opt, 'U'     ,zeros(Mesh_get_numid(mesh), 1));
if ~isfield(opt,'nonlinear') & ~isfield(opt,'kmax'), kmax = 1; end;

if strcmp(nonlinear,'NR')
% -- Conduct Newton-Raphson to find equilibrium position
  Mesh_set_u(mesh,U);
  for k = 1:kmax
    [su,sf,K]  = Mesh_matscale(mesh, Mesh_assemble_k(mesh));
    R  = Mesh_assemble_R(mesh)./sf;
    du = K\R;
    U  = U - du.*su;
    Mesh_set_u(mesh, U);
    fprintf('  %d: %e %e\n', k, norm(du), norm(R));
    if norm(du) < du_tol | norm(R) < R_tol, break; end;
  end
elseif strcmp(nonlinear,'MNR')
% -- Conduct Modified-Newton-Raphson to find equilibrium position
  Mesh_set_u(mesh,U);
  [su,sf,K]  = Mesh_matscale(mesh, Mesh_assemble_k(mesh));
  for k = 1:kmax
    R  = Mesh_assemble_R(mesh)./sf;
    du = K\R;
    U  = U - du.*su;
    Mesh_set_u(mesh, U);
    fprintf('  %d: %e %e\n', k, norm(du), norm(R));
    if norm(du) < du_tol | norm(R) < R_tol, break; end;
  end
else
  error('No such solution method exist');
end

