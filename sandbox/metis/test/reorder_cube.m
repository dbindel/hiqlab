function reorder_cube(ndm,n,prout)
% function reorder_cube(ndm,n,prout)
%
% Reorder matrix from a Laplacian discretized on a cube
%
% Input: ndm   -- dimension of cube
%        n     -- number of nodes per edge
%        prout -- reordering routine
%                 (0:node-based, 1:edge-based)
%

%n  = 15;

e  = ones(n,1);
Ks1= spdiags([-e 2*e -e], -1:1, n, n);
I  = speye(n,n);
if ndm==1
    Ks = Ks1;
elseif ndm==2
    Ks = kron(Ks1,I) + kron(I,Ks1);
elseif ndm==3
    Ks = kron(kron(Ks1,I),I) + kron(kron(I,Ks1),I) + kron(I,kron(I,Ks1));
end

[L,U,P,Q] = lu(Ks);
nnz(L) + nnz(U)

% -- Metis ND
Kst = Ks;
for i = 1:size(Ks,1)
    Kst(i,i) = 0;
end;
[jc,ir] = mat2csc(Kst);
[perm,iperm] = Metis_ND(jc,ir,prout);
perm = perm + ones(size(Kst,1),1);
Knd = Ks(perm,perm);
[Lnd] = chol(Knd);

% -- RCM
[prcm] = symrcm(Ks);
Krcm= Ks(prcm,prcm);
[Lrcm] = chol(Krcm);

% -- AMD
[pamd] = amd(Ks);
Kamd= Ks(pamd,pamd);
[Lamd] = chol(Kamd);

% -- COLAMD
[pcla] = colamd(Ks);
Kcla = Ks(pcla,pcla);
[Lcla] = chol(Kcla);

% -- COLAMD
[Lumf,Uumf,Pumf,Qumf] = lu(Ks);

figure;
subplot(3,2,1);
spy(Lnd+Lnd');
ylabel('Metis')
subplot(3,2,2);
spy(Lrcm+Lrcm');
ylabel('RCM')
subplot(3,2,3);
spy(Lamd+Lamd');
ylabel('AMD')
subplot(3,2,4);
spy(Lcla+Lcla');
ylabel('COLAMD')
subplot(3,2,5);
spy(Lumf+Uumf);
ylabel('UMFPACK')
