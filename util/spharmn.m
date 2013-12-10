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
function Y = spharmn(N, V)
%SPHARMN Creates scalar spherical harmonics of certain degrees.
%
%   Y = SPHARMN(N, V) takes an interval N >= 0 of degrees 
%   represented as a vector and vertices V and computes scalar 
%   spherical harmonics Y_nj of degree n in N and order j=-n,...,n for 
%   every vertex V.
%
%   Note that N must be a vector of consecutive non-negative integers!
%
%   Note Y is of size m-by-g, where m = size(V, 1) is the number of vertices, 
%   g = (N^2 + 2*N - n^2 + 1) is the dimension of the set of scalar 
%   spherical harmonics for degrees N(1),...,N(end).

% Check if N is an interval of consecutive non-negative integers.
assert(isvector(N));
assert(all(N >= 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

m = size(V, 1);

% Compute the dimension.
dim = (N(end)^2 + 2*N(end) - N(1)^2 + 1);

% Create scalar spherical harmonics.
Y = zeros(m, dim);

c = 1;
for k=N
    % Generate scalar spherical harmonics of degree k.
    Y(:, c:(c+2*k)) = spharm(k, V);
    c = c + 2*k + 1;
end

end