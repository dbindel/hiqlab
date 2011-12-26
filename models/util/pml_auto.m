% [beta, N, rmodel] = pml_auto(k, h, order, rtarget)
%
% Compute parameters for a PML to achieve a target
% reflection level.  The PML stretch has the form
%   lambda(x) = 1 - i*beta*x^p
%
% Inputs:
%  k       - vector of wave numbers
%  h       - element size
%  order   - order of element interpolations
%  rtarget - target maximum model reflection (default: 1e-3)
%  p       - polynomial order (default: 2)
%
% Outputs:
%  beta    - coefficient in the PML function
%  N       - number of elements through the PML
%  rmodel  - estimated maximum model reflection

function [beta, N, rmodel] = pml_auto(k, h, order, rtarget, p)

if nargin < 4, rtarget = 1e-3; end
if nargin < 5, p = 2;          end

khmax = max(k)*h;
khmin = min(k)*h;

% -- Find upper and lower bounds on the beta to be considered

beta_l = 1e-7;
beta_u = -(p+1)/2/khmax * log(1e-15);

% -- Bisect to find rinterace = rtol

for j = 1:15
  beta = (beta_u + beta_l)/2;
  ri1 = estimate_pml(khmax, order, beta, p);
  ri2 = estimate_pml(khmin, order, beta, p);
  rinterface = max(ri1, ri2);
  if rinterface < rtarget
    beta_l = beta;
  else
    beta_u = beta;
  end
end

% -- Compute the shortest N that should work

N = ceil(( -(p+1)/2/beta/khmin * log(rtarget) )^( 1/(p+1) ));

% -- Compute the final outputs (rmodel and dimensional beta)

r1 = diagnose_pml1d(khmin, order, beta, p, N);
r2 = diagnose_pml1d(khmax, order, beta, p, N);
rmodel = max(r1, r2);
beta   = beta * h^(-p);

% -- Give a warning if we couldn't get within a factor of 2

if abs(rmodel-rtarget)/rtarget > 1
  warning('Could not get within 2x of target');
end

