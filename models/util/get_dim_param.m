% function [cL,cM,cT,cAux] = get_dim_param(mesh,analysis)
%
%   Function to obtain dimension normalization parameters
%   from a Lua file with the given handle
%Input:
%  mesh      - Mesh object
%  analysis  - 'elastic'       returns [cL,cM,cT]
%            - 'thermoelastic' returns [cL,cM,cT,cTh]
%            - 'electromech'   returns [cL,cM,cT,cA ]
%
function [cL,cM,cT,cAux] = get_dim_param(mesh,analysis)
if nargin < 2, analysis = 'elastic', end;

cL      = Mesh_get_scale(mesh,'L');
cM      = Mesh_get_scale(mesh,'M');
cT      = Mesh_get_scale(mesh,'T');
cAux    = 1;
if strcmp(analysis,'thermoelastic')
   cAux     = Mesh_get_scale(mesh,'Th');
elseif strcmp(analysis,'electromech')
   cAux     = Mesh_get_scale(mesh,'A');
end;
