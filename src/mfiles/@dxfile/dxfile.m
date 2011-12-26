% [dxf] = dxfile(basename)

function [dxf] = dxfile(basename)

dxf.basename = basename;
dxf.fid1     = fopen([basename, '.dx'],  'w');
dxf.fid2     = fopen([basename, '.bin'], 'w', 'b');
dxf.count    = 0;
dxf.objcount = 0;

dxf = class(dxf, 'dxfile');

