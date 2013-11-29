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
function [Y1, Y2] = vspharmn(N, F, V)
%VSPHARMN Creates vector spherical harmonics up to a certain degree.
%
%   [Y1, Y2] = VSPHARMN(N, F, V) takes an interval N > 0 of degrees 
%   represented as a vector and a triangulation F, V and computes vector 
%   spherical harmonics Y_nj^i of degree n in N and order j=-n,...,n for 
%   every face F. Y1, Y2 are the vector spherical harmonics of kind 
%   i={1, 2}.
%
%   Note that N must be a vector of consecutive positive integers!
%
%   Note that the matrices Y1, Y2 contains vectors in R3 and are of size 
%   f-by-g-by-3, where f = size(F, 1) is the number of faces, 
%   g = (N^2 + 2*N - n^2 + 1) is the dimension of the set of vector 
%   spherical harmonics for degrees n,...,N.

% Check if N is an interval of consecutive positive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

n = size(F, 1);

% Compute the dimension of each component.
dim = (N(end)^2 + 2*N(end) - N(1)^2 + 1);

% Create vector spherical harmonics of two kinds.
Y1 = zeros(n, dim, 3);
Y2 = zeros(n, dim, 3);

c = 1;
for k=N
    % Generate vector spherical harmonics of degree k.
    [Yi, Yj] = vspharm(k, F, V);
    Y1(:, c:(c+2*k), :) = Yi;
    Y2(:, c:(c+2*k), :) = Yj;
    c = c + 2*k + 1;
end

end