function gen_gauss_tester(fname)

  fp = fopen(fname, 'w');
  fprintf(fp, '#include <cstdio>\n');
  fprintf(fp, '#include <cmath>\n');
  fprintf(fp, '#include \"gaussquad.h\"\n');
  fprintf(fp, 'int main()\n');
  fprintf(fp, '{\n');
  for N = 1:10
    [x,w] = gauss(N);
    for j = 1:N
      fprintf(fp, '    if (fabs(gauss_point(%d,%d)-(%0.16f)) > 1e-15)\n', ...
              j-1,N,x(j));
      fprintf(fp, '        printf(\"Failure on abscissa (%d,%d)\\n\");\n',j,N);
      fprintf(fp, '    if (fabs(gauss_weight(%d,%d)-(%0.16f)) > 1e-15)\n', ...
              j-1,N,w(j));
      fprintf(fp, '        printf(\"Failure on abscissa (%d,%d)\\n\");\n',j,N);
    end
  end
  fprintf(fp, '}\n');


% (From LNT's /Spectral Methods in MATLAB/)
% GAUSS  nodes x (Legendre points) and weights w
%        for Gauss quadrature

function [x,w] = gauss(N)

  beta = .5./sqrt(1-(2*(1:N-1)).^(-2));
  T = diag(beta,1) + diag(beta,-1);
  [V,D] = eig(T);
  x = diag(D); [x,i] = sort(x);
  w = 2*V(1,i).^2;
