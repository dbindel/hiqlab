function y = spdiag(x)

if isvector(x)
  y = spdiags(x,0,length(x),length(x));
else
  y = full(diag(x));
end
