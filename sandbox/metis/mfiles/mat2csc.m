function [jc,ir,val] = mat2csc(A)
% function [jc,ir,val] = mat2csc(A)
%
% Converts Matrix to CSC format
%
[ir,c,val] = find(A);

n   = size(A,1);
nnz = length(ir);
jc  = zeros(n+1,1);

% -- Assume results in CSC format
jc(1) = 1;
for i=1:n
   ind  = find(c==i);
   jc(i+1)= jc(i) + length(ind);
end;

% -- Convert to 0-base
jc = jc - ones(n+1,1);
for i = 1:nnz
    ir(i) = ir(i) - 1;
end
