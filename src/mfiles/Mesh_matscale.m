function [su,sf,varargout] = Mesh_matscale(mesh,varargin)

[su,sf] = Mesh_get_scaling_vectors(mesh);
n = length(su);
Du = spdiag(su);
Df = spdiag(1./sf);

nout = nargin-1;
for i = 1:nout
  varargout{i} = Df*varargin{i}*Du;
end
