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
function d = vspharmeigs(N)
%VSPHARMEIGS Computes a vector of eigenvalues for vector spherical
%harmonics.
%
%   d = VSPHARMEIGS(N) takes degrees N and returns a vector of eigenvalues
%   of the vector spherical harmonics.
%
%   Note that N must be an interval of consecutive positive integers!
%
%   d is a vector of length dim, where dim is the dimension of the vector
%   spherical harmonics induced by degrees N.

% Check if N is an interval of positive consecutive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

% Compute eigenvalues.
d = spharmeigs(N);

% Duplicate since eigenvalues are same for first and second kind.
d = repmat(d, 2, 1);

end