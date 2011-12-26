% x = qqoptdefault(opt, name, default)
% Get the named parameter value from the opt structure.  If the parameter
% doesn't exist in the structure, set it to the default value.

function x = qoptdefault(opt, name, default)

% HiQLab
% Copyright (c): Regents of the University of California

if isfield(opt, name)
  x = getfield(opt, name);
else
  x = default;
end
