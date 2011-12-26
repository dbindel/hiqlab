 
% Test case: Q convergence sweep

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_converge2.m,v 1.1 2006/05/01 05:04:09 tkoyama Exp $
if 0
wall = [];
Qall = [];
tall = [];

for order = 1:3
  for mdense = 1:8

    tic; 
    param.order = order;
    param.dense = 1e-6/mdense;
    [w,Q] = driver_mode(2*pi * 715e6, param);
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
end
set(0, 'defaulttextinterpreter', 'none');

figure(1);
iall = 1:8;
h = plot(iall, Qall(1,:), '-', ...
     iall, Qall(2,:), ':', ...
     iall, Qall(3,:), '-.');
set(h, 'LineWidth', 2);
legend('Linear', 'Quadratic', 'Cubic', 0);
ylabel('Computed $Q$');
xlabel('Mesh density');
laprint(1, 'converge_Q', 'width', 12);

figure(2);
iall = 1:8;
h = plot(iall, real(wall(1,:))/1e6, '-', ...
     iall, real(wall(2,:))/1e6, ':', ...
     iall, real(wall(3,:))/1e6, '-.');
set(h, 'LineWidth', 2);
legend('Linear', 'Quadratic', 'Cubic', 0);
ylabel('Computed frequency (MHz)');
xlabel('Mesh density');
laprint(2, 'converge_w', 'width', 12);
