
function [dxf] = writefield(dxf, varargin)

dxf.objcount = dxf.objcount + 1;

if mod(length(varargin),2) == 0
  fprintf(dxf.fid1, 'object %d class field\n', dxf.objcount);
else
  fprintf(dxf.fid1, 'object "%s" class field\n', varargin{1});
  varargin = varargin(2:end);
end

for k = 1:2:length(varargin)
  fprintf(dxf.fid1, 'component "%s" value %d\n', ...
	  varargin{k}, varargin{k+1});
end
fprintf(dxf.fid1, '#\n');
