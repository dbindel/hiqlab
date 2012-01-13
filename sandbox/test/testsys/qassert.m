function qassert(flag, msg)

if ~flag 
  error(sprintf('Failed: %s', msg));
end
