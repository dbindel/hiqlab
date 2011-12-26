 
% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_all.m,v 1.1 2006/05/01 05:04:09 tkoyama Exp $
 
% -- By default run everything
if ~exist('which_test')
  which_test = 0;
end

% -- Run a convergence test
if (which_test == 0) | (which_test == 1)
  test_converge;
  print('-deps', 'Qconverge.eps');
end

% -- Run the Bode plot comparison
if (which_test == 0) | (which_test == 2)
  test_rombode;
  print('-deps', 'rombode1.eps');
end

% -- Simple eigensweep: Q and root locus
if (which_test == 0) | (which_test == 3)
  test_sweep0;
  print('-f1', '-deps', 'sweeptest0a.eps');
  print('-f2', '-deps', 'sweeptest0b.eps');
  print('-f3', '-deps', 'sweeptest0c.eps');
end

% -- Compare Q vs simulation for three PML settings
if (which_test == 0) | (which_test == 3)
  test_sweep1;
  print('-deps', 'sweeptest1.eps');
end

% -- Movie sweep
if (which_test == 0) | (which_test == 3)
  test_sweep2;
end

% -- Diamond data sweep
if (which_test == 0) | (which_test == 4)
  test_diamond;
end
