function partition_cube(ndm,n,nparts,prout)
% function partition_cube(ndm,n,nparts,prout)
% 
% Partition a ndm-cube into nparts with Metis
%
% Input: ndm    -- dimensions of cube
%        n      -- nodes per edge
%        nparts -- number of parts to split
%        prout  -- partitioning routine
%                  (0: recursive bisection, 1:k-way)
%
if nargin < 4, prout = 0; end

% -- Parameters
%ndm   = 3;   % -- Dimension of cube
%n     = 10;  % -- Nodes per side of the cube
%nparts= 4;   % -- Number of partitions
%prout = 0;   % -- Partitioning routing (0:Recursive bisection, 1:K-Way)

% -- Construct matrix representation
e  = ones(n,1);
Ks1= spdiags([-e 2*e -e], -1:1, n, n);
I  = speye(n,n);
x1d= linspace(0,1,n);
e1d= ones(1,n);

if ndm==1
    Ks = Ks1;
    x  = x1d;
elseif ndm==2
    Ks = kron(Ks1,I) + kron(I,Ks1);
    x  = kron(x1d,e1d);
    y  = kron(e1d,x1d);
elseif ndm==3
    Ks = kron(kron(Ks1,I),I) + kron(kron(I,Ks1),I) + kron(I,kron(I,Ks1));
    x  = kron(kron(x1d,e1d),e1d);
    y  = kron(kron(e1d,x1d),e1d);
    z  = kron(kron(e1d,e1d),x1d);
end

% -- Construct xadj and adj
for i = 1:size(Ks,1)
    Ks(i,i) = 0;
end;
[jc,ir] = mat2csc(Ks);

% -- Partition with METIS
[part,edgecut] = Metis_PartGraph(nparts,jc,ir,prout);

% -- Plot
strc = {'r','b','g','m','c','y','k'};
strm = {'*','x','o','+','s','d','p'};

figure;
if ndm==1

    e1d= zeros(1,n);
    hold on;
    for i = 1:nparts
        ids = find(part==(i-1));
        ic  = mod(i-1,length(strc))+1;
        im  = mod(floor(i/length(strc)),length(strm))+1;
        plot(x(ids),e1d(ids),strcat(strc{ic},strm{im}));
    end
    hold off;
    axis equal;
    axis([-0.5,1.5,-0.5,0.5]);

elseif ndm==2

    hold on;
    for i = 1:nparts
        ids = find(part==(i-1));
        ic  = mod(i-1,length(strc))+1;
        im  = mod(floor(i/length(strc)),length(strm))+1;
        plot(x(ids),y(ids),strcat(strc{ic},strm{im}));
    end
    hold off;
    axis equal;
    axis([-0.5,1.5,-0.5,1.5]);

elseif ndm==3

    hold on;
    for i = 1:nparts
        ids = find(part==(i-1));
        ic  = mod(i-1,length(strc))+1;
        im  = mod(floor(i/length(strc)),length(strm))+1;
        plot3(x(ids),y(ids),z(ids),strcat(strc{ic},strm{im}));
    end
    hold off;
    axis equal;
    axis([-0.5,1.5,-0.5,1.5,-0.5,1.5]);

end
