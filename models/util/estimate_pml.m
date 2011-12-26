
% [rinterface, rcontinuum] = estimate_pml(k, order, c, p, N)
%
% Estimate the interface reflection and the continuum reflection
% given a set of PML parameters.
%
% Inputs:
%   k      - desired wave number
%   order  - element order
%   c      - PML profile coefficient
%   p      - PML profile exponent
%   N      - number of elements in the PML
%
% Outputs:
%   rinterface - interface reflection coefficient
%   rcontinuum - continuum reflection coefficient

function [rinterface, rcontinuum] = estimate_pml(kx, order, c, p, N)

N2 = ceil( (-(p+1)/2/c/kx * log(eps))^(1/(p+1)) );
if nargin < 5, N = N2; end

rinterface = diagnose_pml1d(kx, order, c, p, N2);
rcontinuum = exp(-2*c*kx/(p+1) * N^(p+1));

