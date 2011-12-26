% [V,w,Q] = emmode(mesh, w0, nev, opt);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
%
% opt   contains optional parameters:
%   mech - Conduct purely mechanical analysis
%   use_matlab - use matlab eigs or not (Default:1)
%   idg_m - array of global numbers for mechanical variables
%   idg_p - array of global numbers for potential variables

function [V,w,Q] = emmode(mesh, w0, nev, opt);

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 3, nev = 1;    end
if nargin < 4, opt = []; end

mech           = qoptdefault(opt, 'mech',      0);
opt.use_matlab = qoptdefault(opt, 'use_matlab', 1);
opt.idg_m      = qoptdefault(opt, 'idg_m', []);
opt.idg_p      = qoptdefault(opt, 'idg_p', []);

if mech
   [V,w,Q] = emmode_mech_1(mesh,w0,nev,opt);
   return
end

if opt.use_matlab
  [V,w,Q] = emmode_matlab_1(mesh, w0, nev, opt);
end


% -------------------------------------------------------------------
% [V,w,Q] = emmode_matlab_1(mesh, w0, nev, opt);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
% (mechanical and thermal parts).
%
% opt contains optional parameters:
%   idg_m - array of global numbers for mechanical variables
%   idg_p - array of global numbers for potential variables
%
function [V,w,Q] = emmode_matlab_1(mesh, w0, nev, opt);

if nargin < 3, nev = 1;    end
if nargin < 4, opt = []; end

idg_m = qoptdefault(opt, 'idg_m', []);
idg_p = qoptdefault(opt, 'idg_p', []);

numid   = Mesh_get_numid(mesh);    % Number of IDs

% -- Extract nodal ids
numid  = Mesh_get_numid(mesh);
id     = Mesh_get_id(mesh);
ndm    = Mesh_get_ndm(mesh);
ndf    = Mesh_get_ndf(mesh);
id_m   = [];

% -- Number of nodal mechanical dofs
for i = 1:ndf-1
    tmp = find(id(i,:)>0);
    id_m= union(id_m,id(i,tmp));
end
% -- Number of nodal potential dofs
id_p  = find(id(ndf,:)>0);
id_p  = id(ndf,id_p);

% -- Extract global variable ids
for i = 1:length(idg_m)
    idg_m(i) = Mesh_globalid(mesh,idg_m(i));
end
for i = 1:length(idg_p)
    idg_p(i) = Mesh_globalid(mesh,idg_p(i));
end

% -- Sort ids
id_d1=[id_m,idg_m,id_p,idg_p];

nid_d1 = length(id_d1);

% -- Extract matrices
[M,K,su] = Mesh_assemble_mk_nd(mesh);
[V,w,Q]  = pml_mode(M(id_d1,id_d1),K(id_d1,id_d1),w0,nev,opt);
V(id_d1,:) = V;
V = spdiag(su)*V;

V = V / diag(sqrt(diag(V'*V)));
Q = abs(w)./(2*imag(w));


% -------------------------------------------------------------------
% [V,w,Q] = emmode_mech_1(mesh, w0, nev, opt);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
% (thermal part is zero) Compute only the mechanical frequency

function [V,w,Q] = emmode_mech_1(mesh, w0, nev, opt);

if nargin < 3, nev = 1;    end
if nargin < 4, opt = []; end

idg_m= qoptdefault(opt, 'idg_m', []);
idg_p= qoptdefault(opt, 'idg_p', []);

numid    = Mesh_get_numid(mesh);      % Number of IDs
id       = Mesh_get_id(mesh);
ndf      = Mesh_get_ndf(mesh);
id_m     = [];

% -- Number of nodal mechanical dofs
for i = 1:ndf-1
    tmp = find(id(i,:)>0);
    id_m= union(id_m,id(i,tmp));
end

% -- Number of global mechanical dofs
for i = 1:length(idg_m)
    idg_m(i) = Mesh_globalid(mesh,idg_m(i));
end
id_m        = union(id_m,idg_m);
nm           = length(id_m);
[Muu,Kuu,su] = Mesh_assemble_mk_nd(mesh);
Muu          = Muu(id_m,id_m);
Kuu          = Kuu(id_m,id_m);

[VV,w,Q]  = pml_mode(Muu,Kuu,w0,nev);
V         = zeros(numid,nev);
V(id_m,:) = VV;
V         = spdiag(su) * V;
