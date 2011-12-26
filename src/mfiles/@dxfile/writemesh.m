% [dxf] = writemesh(dxf, p, e)

function [dxf] = writemesh(dxf, p, e)

if size(p,2) == 2
  dxf = writemesh2d(dxf, p, e);
elseif size(p,2) == 3
  dxf = writemesh3d(dxf, p, e);  
else
  error('Wrong size');
end


function [dxf] = writemesh2d(dxf, p, e)

eplot = plotelt2d(e);
dxf = writearray(dxf, p, 'float', ...
		 'dep', 'positions');
dxf = writearray(dxf, eplot([1 2 4 3],:)-1, 'int', ...
		 'element type', 'quads',    ...
		 'ref',          'positions' );


function [dxf] = writemesh3d(dxf, p, e)

if size(e,1) == 8,       m = 1;
elseif size(e,1) == 27,  m = 2;
elseif size(e,1) == 64,  m = 3;
else,                    error('Bad element array size');
end

eplot = [];
dx = (m+1)*(m+1);
dy = (m+1);
dz = 1;

for ix = 0:m-1
  for iy = 0:m-1
    for iz = 0:m-1
      
      eplot = [eplot, ...
	       e([(ix  )*dx + (iy  )*dy + (iz  ), ...
		  (ix  )*dx + (iy  )*dy + (iz+1), ...
		  (ix  )*dx + (iy+1)*dy + (iz  ), ...
		  (ix  )*dx + (iy+1)*dy + (iz+1), ...
		  (ix+1)*dx + (iy  )*dy + (iz  ), ...
		  (ix+1)*dx + (iy  )*dy + (iz+1), ...
		  (ix+1)*dx + (iy+1)*dy + (iz  ), ...
 		  (ix+1)*dx + (iy+1)*dy + (iz+1)] + 1, :)];
    end
  end
end

dxf = writearray(dxf, p, 'float', ...
		 'dep', 'positions');
dxf = writearray(dxf, eplot-1, 'int', ...
		 'element type', 'cubes',    ...
		 'ref',          'positions' );
