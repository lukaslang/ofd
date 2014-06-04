% Copyright 2013 Clemens Kirisits and Lukas Lang
%
% This file is part of OFD.
%
%    OFD is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    OFD is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with OFD.  If not, see <http://www.gnu.org/licenses/>.
function [S, rho] = surfsynth(N, X, c, Y)
%SURFSYNTH Computes synthesis of a sphere-like surface for points on the
%sphere.
%
%   [S, rho] = surfsynth(N, X, c) takes coefficients c for scalar sphercial 
%   harmonics of degrees N and points X on the unit sphere and returns the
%   projection V of these to the sphere-like surface. rho is the radial
%   function evaluated at X.
%
%   [S, rho] = surfsynth(N, X, c, Y) additionally takes scalar spherical harmonics
%   Y which have been evaluated at X. This is useful e.g. after calling 
%   surffit, which returns Y.
%
%   Note that N must be a vector of non-negative consecutive integers. 
%   X is an n-by-3 matrix of points on the unit sphere. c is a vector of 
%   size dim, where dim is the number of scalar spherical harmonics of 
%   degrees specified by N.
%
%   S is a matrix of size n-by-3. rho is a vector of length n.

% Check if N is an interval of consecutive positive integers.
assert(isvector(N));
assert(all(N >= 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

% Compute and check dimension.
dim = N(end)^2 + 2*N(end) - N(1)^2 + 1;
assert(isvector(c));
assert(length(c) == dim);

if(nargin == 4)
    % Check if dimensions comply.
    assert(size(Y, 1) == size(X, 1));
    assert(size(Y, 2) == dim);
else
    % Evaluate scalar spherical harmonics at points X.
    Y = spharmn(N, X);
end

% Recover surface function at X.
rho = Y * c;

% Compute coordinates of surface at points X.
S = bsxfun(@times, X, rho);