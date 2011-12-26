% [V,w,Q] = tedmode(mesh, w0, nev, opt);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
% (mechanical and thermal parts).
%
% *Note*: tedmode will reorder the degrees of freedom in the mesh.
%
% opt   contains optional parameters:
%   mech - Conduct purely mechanical analysis
%   type - Use a perturbation method ('pert') or linearization
%          ('full').  Default is 'pert'.
%   T0   - Baseline temperature for use in symmetric linearization
%   use_matlab - use matlab eigs or not (Default:0)

function [V,w,Q] = tedmode(mesh, w0, nev, opt);

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 3, nev = 1;  end
if nargin < 4, opt = []; end

mech           = qoptdefault(opt, 'mech',       0     );
typ            = qoptdefault(opt, 'type',       'pert');
T0             = qoptdefault(opt, 'T0',         0     );
opt.use_matlab = qoptdefault(opt, 'use_matlab', 0     );

if mech
  [V,w,Q] = tedmode_mech_1(mesh,w0,nev,opt);
  return
end

if opt.use_matlab
  [V,w,Q] = tedmode_matlab_1(mesh, w0, nev, opt);
else
  cT      = Mesh_get_scale(mesh, 'T'); % Characteristic time scale
  numid   = Mesh_get_numid(mesh);    % Number of IDs
  nm      = ted_block_mesh(mesh);    % Number of mech IDs
  nt      = numid - nm;              % Number of thermal IDs
  ncv     = max(10,2*nev);
  if strcmp(typ, 'pert')
    [w,V] = ted_compute_eigsp(mesh, w0, nev, ncv);
  else %if strcmp(typ, 'full')
    [w,V] = ted_compute_eigs(mesh, w0, T0, cT, nev, ncv);
  end;
  V = V / diag(sqrt(diag(V'*V)));
  Q = abs(w)./(2*imag(w));
end


% ------------------------------------------------------------------
% [V,w,Q] = tedmode_matlab_1(mesh, w0, nev, param);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
% (mechanical and thermal parts).
%
% *Note*: tedmode will reorder the degrees of freedom in the mesh.
%
% param contains optional parameters:
%   type - Use a perturbation method ('pert') or linearization
%          ('full').  Default is 'pert'.
%   T0   - Baseline temperature for use in symmetric linearization

function [V,w,Q] = tedmode_matlab_1(mesh, w0, nev, param);

if nargin < 3, nev = 1;    end
if nargin < 4, param = []; end

typ = qoptdefault(param, 'type', 'pert');
T0  = qoptdefault(param, 'T0',   0     );

numid = Mesh_get_numid(mesh);    % Number of IDs
nm    = ted_block_mesh(mesh);    % Number of mech IDs
nt    = numid - nm;              % Number of thermal IDs
[M,K,C, su] = Mesh_assemble_mkc_nd(mesh);

Iu = 1:nm;
It = nm+1:numid;

% -- Extract system submatrices
Kuu = K(Iu,Iu);  % Mechanical stiffness
Muu = M(Iu,Iu);  % Mechanical mass
Ctt = C(It,It);  % Thermal mass
Ktt = K(It,It);  % Thermal diffusion term
Ctu = C(It,Iu);  % Mechanically driven heating
Kut = K(Iu,It);  % Thermally driven expansion

% -- Symmetrize matrices that should be symmetric
Kuu = (Kuu + Kuu.')/2;
Muu = (Muu + Muu.')/2;
Ctt = (Ctt + Ctt.')/2;
Ktt = (Ktt + Ktt.')/2;

if strcmp(typ, 'zener')

  % -- Preliminary eigencomputation on undamped system
  pmlopt = param;
  pmlopt.use_matlab = 1;
  [U,D] = pml_mode(Muu,Kuu,w0,nev,pmlopt);
  w0    = D;

  for k = 1:length(w0)
    % -- Compute the resulting temperature fields
    Th(:,k) = -(Ctt - complex(0,1)/w0(k)*Ktt)\(Ctu*U(:,k));

    % -- Solve an update equation
    dw2  = ( U(:,k)'*(Kut*Th(:,k)) )/( U(:,k)'*Muu*U(:,k) );
    w(k) = sqrt(w0(k)^2 + dw2);
    Q(k) = abs(w(k))/(2*imag(w(k)));

  end

  V = [U; Th];
  V = spdiag(su) * V;
  w = w.';
  Q = Q.';

elseif strcmp(typ, 'pert')

  % -- Preliminary eigencomputation on undamped system
  pmlopt.use_matlab = 1;
  [U,D] = pml_mode(Muu,Kuu,w0,nev,pmlopt);
  w0    = D;

  for k = 1:length(w0)
    % -- Compute the resulting temperature fields
    Th(:,k) = -(Ctt - complex(0,1)/w0(k)*Ktt)\(Ctu*U(:,k));

    % -- Solve bordered system for perturbation
    delta = [ Kuu - w0(k)^2*Muu, -2*w0(k)*Muu*U(:,k);
              1e-7*U(:,k)',    0              ] \ [-Kut*Th(:,k); 0];
    du    = delta(1:end-1);
    dw    = delta(end);

    % -- Apply the perturbation
    U(:,k) = U(:,k)+du;
    w(k)   = w0(k)+dw;
    Q(k)   = abs(w(k))/(2*imag(w(k)));

  end

  V = [U; Th];
  V = spdiag(su) * V;
  w = w.';
  Q = Q.';

else %if strcmp(typ, 'full')

  cT = Mesh_get_scale(mesh, 'T');
  Zuu  = spalloc(nm,nm,0);
  Zut  = spalloc(nm,nt,0);
  Iuu  = speye(nm);

  B    = [  Iuu ,  Zuu ,       Zut ;
            Zuu ,  Muu/cT^2 ,  Zut ;
            Zut',  Zut',       Ctt/cT ];

  A    = [  Zuu ,  Iuu ,     Zut ;
           -Kuu ,  Zuu ,    -Kut ;
            Zut', -Ctu/cT , -Ktt ];

  pmlopt = param; % DSB
  pmlopt.linearized = 1;
  pmlopt.use_matlab = 1;
  [V,w,Q] = pml_mode(B,A,1i*w0*cT,nev,pmlopt);
  w = w/1i/cT;
  V = V([1:nm, 2*nm+1:2*nm+nt],:);
  V = spdiag(su) * V;
  Q = abs(w)./(2*imag(w));

end


% ------------------------------------------------------------------
% [V,w,Q] = tedmode_mech_1(mesh, w0, nev, opt);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
% (thermal part is zero) Compute only the mechanical frequency
%
% *Note*: tedmode will reorder the degrees of freedom in the mesh.

function [V,w,Q] = tedmode_mech_1(mesh, w0, nev, opt);

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 3, nev = 1;    end
if nargin < 4, opt = []; end

numid = Mesh_get_numid(mesh);      % Number of IDs
nm    = ted_block_mesh(mesh);      % Number of mech IDs

[Muu,Kuu, su] = Mesh_assemble_mk_nd(mesh);
Muu = Muu(1:nm,1:nm);
Kuu = Kuu(1:nm,1:nm);

[VV,w,Q]  = pml_mode(Muu,Kuu,w0,nev);
V         = zeros(numid,nev);
V(1:nm,:) = VV;
V         = spdiag(su(1:nm))*VV;

