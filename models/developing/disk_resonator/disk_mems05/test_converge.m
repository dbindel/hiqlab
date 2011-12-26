 
% Test case: Q convergence sweep

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_converge.m,v 1.1 2006/05/01 05:06:55 tkoyama Exp $

wall = [];
Qall = [];
tall = [];

for order = 1:3
  for mdense = 1:4

    tic; 
    param.order = order;
    param.dense = 1e-6/mdense;
    [w,Q] = driver_mode(2*pi * 47.5e6, param);
    w = w/2/pi;
    t = toc;

    fprintf('\r%d (%d):\t%4.1f MHz\t%4.0f Q\t%4.1f s', ...
             order, mdense, real(w)/1e6, Q, t);
    wall(order, mdense) = w;
    Qall(order, mdense) = Q;
    tall(order, mdense) = t;

  end
end
fprintf('\n');

plot(Qall');
legend('Linear', 'Quadratic', 'Cubic', 0);
ylabel('Computed Q');
xlabel('Mesh density');
