% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_ffbeam3d.m,v 1.1 2006/05/01 07:45:36 tkoyama Exp $
%
% [freq ,Q, mesh, Lua_handle] = driver_ffbeam3d(params)
%
%   Computes Frequency and Q factor
%
%   Input: *params        -Additional analysis parameters
%           .analysisType -'full' for Full Eigen computation
%                          'pert' for Perturbation computation
%           .use_matlab   -'1' to use matlab solvers
%                         -'0' to ship outside of matlab
%           .order        -Order of interpolation of elements
%           .dense        -Density of mesh
%           .f0           -Stretch parameter for PML
%           .freqN        -Shift for eigen computation
%   Remark: '*' are optional parameters
%
%

function [freq ,Q, mesh, Lua_handle] = driver_ffbeam3d(params)
 
if nargin < 1, params = [];      end;

numeigs   = 1;

use_matlab   = getparam(params, 'use_matlab', 0);
params.order = getparam(params, 'order',  2);
params.dense = getparam(params, 'dense', 10);
analysisType = getparam(params, 'analysisType', 'full');
f0           = getparam(params, 'f0'    ,    20); 
freqN        = getparam(params, 'shift' , 9.6e6);
w0           = freqN*2*pi;

disp('Loading Lua mesh file');
meshfile            = 'ffbeam3d_mesh.lua';
[mesh, Lua_handle]  = Mesh_load(meshfile, params);
[cL,cM,cT,cTh]      = get_dim_param (Lua_handle,'thermoelastic');
T0                  = Lua_get_double(Lua_handle,'T0');
numid               = Mesh_get_numid(mesh);
fprintf('Number of id: %g\n',numid);
disp('Go Ahead?');
pause;


disp('** Thermoelastic-PML **');
p_full.cT         = cT;
p_full.T0         = T0;
p_full.type       = analysisType;
p_full.use_matlab = use_matlab;
[V,w,Q]     = tedmode(mesh,w0,numeigs,p_full);
freq        = w/2/pi;

% -- Print results
fprintf('Original frequencey    : %g\n',freqN);
fprintf('Damped  frequency(real): %g\n',real(freq));
fprintf('Damped  frequency(imag): %g\n',imag(freq));
fprintf('Shift(w0)              : %g\n',w0/2/pi);
fprintf('Damped  w        (real): %g\n',real(w)/2/pi);
fprintf('Damped  w        (imag): %g\n',imag(w)/2/pi);
fprintf('Q                      : %g\n',Q);
fprintf('Numid                  : %g\n',numid);

% -- Put back mode vector and plot mesh
Mesh_set_u(mesh, V(:,1));

% -- Plot mesh
ndof = Mesh_get_ndf(mesh);
if isfield(params,'plotmesh')
  if params.plotmesh
    figure(1);
    if ndof == 3 
      plotmesh2d(meshfile,params);
      axis equal;
    elseif ndof == 4
      plotmesh3d(meshfile,params);
    end;  
  end;
end;

% -- Plot deformed shape
if isfield(params,'plotdefo')
  if params.plotdefo  
    figure(2);
    plotopt = [];
    if ndof == 3
      plotopt.cfields = [3];
      plotopt.axequal = 1;
      plotopt.deform  = getparam(params , 'deform' ,  5                  );
      plotopt.deformc = getparam(params , 'deformc',  5                  );
      plotopt.n_div   = getparam(params , 'n_div'  , 20);
      plotopt.n_iter  = getparam(params , 'n_iter' , 60);
      plotfield2d(mesh,Lua_handle,plotopt);
    elseif ndof == 4 
      s_axis          = getparam(params , 's_axis' ,  3   );
      s_coord         = getparam(params , 's_coord', 'z1' );
      plotopt.s_axis  = s_axis;
      plotopt.s_coord = s_coord;
      plotopt.s_coord = Lua_get_double(Lua_handle,s_coord)*cL;

      plotopt.ufields = getparam(params , 'ufields', setdiff(1:3,s_axis) );
      plotopt.cfields = getparam(params , 'cfields', [4]                 );
      plotopt.axequal = 1;
      plotopt.deform  = getparam(params , 'deform' ,  5                  );
      plotopt.deformc = getparam(params , 'deformc',  5                  );
      plotopt.n_div   = getparam(params , 'n_div'  , 20);
      plotopt.n_iter  = getparam(params , 'n_iter' , 60);
      plotfield3d(mesh,Lua_handle,plotopt);
    end;
  end;
end;

if nargout < 3
  Mesh_delete(mesh);
end;


% ---------------------------------------------------------------
%
% Get named parameter or default value.
%
function [v] = getparam(p, name, default)

if isfield(p, name)
    v = getfield(p, name);
else
    v = default;
end
