function [A] = csc2mat(xadj,adj)
% function [A] = csc2mat(xadj,adj)
% 
% Converts CSC to Matrix format
%
nv = length(xadj) - 1;
nnz= length(adj);
A  = spalloc(nv,nv,nnz);

for i = 1:nv

    for j = 1:(xadj(i+1)-xadj(i))

        A(adj(xadj(i)+j)+1,i) = 1;

    end

end
