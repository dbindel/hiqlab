% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_piezo_capacitor2d_f_bode.m,v 1.1 2006/05/01 08:19:07 tkoyama Exp $
%
% function [H,freq] = drive_em_capacitor2d_f_bode(meshfile,params,w0,opt)
%
%   Plots forced response
%
%   Input:  meshfile      -Filename of the mesh input
%          *params        -Additional analysis parameters
%           w0            - center frequency[Hz]
%          *opt
%          -plot_mesh(0)  -'1' to plot mesh
%          -wr_min(0.90)  - left  value mag of bode plot
%          -wr_max(1.10)  - right value mag of bode plot
%          -w_ndiv(50)    - number divisions in bode plot
%          -kmax(0)       - number of arnoldi iterations
%
%
function [H,freq] = driver_em_capacitor2d_f_bode(meshfile,params,w0,opt)

if nargin < 2, error('Need Lua mesh input file'); end;
if nargin < 3, w0 = params;     params = [];      end;
if nargin < 4
  if ~isnumeric(w0), error('Need center driving frequency'); 
  else,            opt    = [];                     end;
end

f0                 = qoptdefault(params, 'f0'       , 0.0);

plot_mesh = qoptdefault(opt, 'plot_mesh',    0);
wr_max    = qoptdefault(opt, 'wr_max'   , 1.10);
wr_min    = qoptdefault(opt, 'wr_min'   , 0.90);
w_ndiv    = qoptdefault(opt, 'w_ndiv'   ,   50);

[mesh, Lua_handle] = Mesh_load(meshfile,params);
[cL, cM, cT, cA  ] = get_dim_param( mesh, 'electromech');
cR                 = Mesh_get_scale(mesh, 'R');
numid              = Mesh_get_numid(mesh);
nm                 = em_block_mesh(mesh);
wcenter            = w0*2*pi*cT;
freq               = linspace(wcenter*wr_min,wcenter*wr_max,w_ndiv);

fprintf('#Free DOF      : %g\n',numid);
fprintf('#Free DOF(mech): %g\n',nm);
fprintf('Proceed?\n');
plotmesh(mesh);
axis equal;
pause;

% -- Get info on Top and Bottom(TieField) elements
% -- Top
Vt_a    = Lua_get_double(Lua_handle,'Vt_a');
ntie_t  = Lua_get_double(Lua_handle,'ntie_t')+1;
ntie_tI = Mesh_branchid(mesh,1,ntie_t);
ntie_tA = Mesh_nbranch_idja(mesh,ntie_t);
ntie_tE = Mesh_branchid(mesh,ntie_tA,ntie_t);
% -- Bottom
Vb_a    = Lua_get_double(Lua_handle,'Vb_a');
ntie_b  = Lua_get_double(Lua_handle,'ntie_b')+1;
ntie_bI = Mesh_branchid(mesh,1,ntie_b);
ntie_bA = Mesh_nbranch_idja(mesh,ntie_b);
ntie_bE = Mesh_branchid(mesh,ntie_bA,ntie_b);

% -- Solve for AC fluctuation
% -- Relabel ID's and form matrices
[M,K]              = Mesh_assemble_mk(mesh);
subplot(1,2,1);
spy(M);
subplot(1,2,2);
spy(K);

F                  = zeros(numid,1);
F(ntie_tE,1)       = Vt_a;
F(ntie_bE,1)       = Vb_a;

% -- Compute transfer function
kmax = 0;
use_umfpack = 1;
use_umfpack = (use_umfpack & exist('umfpack'));

H = zeros(length(freq),1);

for idx = 1:length(freq)
  w = freq(idx);
  if (kmax>=0)
    fprintf('%d: %d Hz\n', idx, freq(idx)/cT/2/pi);
  end;

  if use_umfpack
    U = umfpack(K - w^2*M, '\', F);
  else          
    U = (K - w^2*M)\F;
  end
  Vt= U(ntie_tI,1);
  Qt= U(ntie_tE,1);
  H(idx) = Qt*complex(0,w)/Vt*cA;
%  H(idx) = Qt*complex(0,w)/Vt/cR;
end

% -- Plot mesh
ndof = Mesh_get_ndf(mesh);
if plot_mesh
  figure(1);
  if ndof == 3
    plotmesh(mesh);
    axis equal;
  elseif ndof == 4
    plotmesh3d(meshfile,params);
  end;
end;
Mesh_delete(mesh);

% -- Bode Plot
fprintf('\n');
fprintf('------------------Bode Plot----------------\n');
fprintf('PML f0        : %g\n',f0);
fprintf('Center F  [Hz]: %g\n',w0);
fprintf('No. of id     : %g\n',numid);

figure(2)
freq = freq/cT/2/pi;
opt.lstyle = 'r*-';
opt.visualQ= 1;
plot_bode(freq,H,opt);
