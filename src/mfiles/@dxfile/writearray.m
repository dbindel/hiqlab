function [dxf] = writearray(dxf, A, type, varargin)

fid1 = dxf.fid1;
fid2 = dxf.fid2;
dxf.objcount = dxf.objcount + 1;

if (strcmp(type, 'float'))
  filetype = 'real*4';
elseif (strcmp(type, 'int'))
  filetype = 'integer*4';
else
  error('Invalid dx array type');
end

fprintf(fid1, 'object %d class array type %s\n', dxf.objcount, type);

if size(A,1) == 1 | size(A,2) == 1
  fprintf(fid1, 'rank 0 items %d\n', length(A));
else
  fprintf(fid1, 'rank 1 shape %d items %d\n', size(A,1), size(A,2));
end

fprintf(fid1, 'msb binary data file %s.bin,%d\n', ...
	dxf.basename, dxf.count*4);

for k = 1:2:length(varargin)
  fprintf(fid1, 'attribute "%s" string "%s"\n', varargin{k}, varargin{k+1});
end

fprintf(fid1, '#\n');

dxf.count = dxf.count + fwrite(fid2, A, filetype);

