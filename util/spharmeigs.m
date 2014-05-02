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
function d = spharmeigs(N)
%SPHARMEIGS Computes a vector of eigenvalues for scalar spherical
%harmonics.
%
%   d = SPHARMEIGS(N) takes degrees N and returns a vector of eigenvalues
%   of scalar spherical harmonics.
%
%   Note that N must be an interval of consecutive non-negative integers!
%
%   d is a vector of length dim, where dim is the dimension of the scalar
%   spherical harmonics induced by degrees N.

% Check if N is an interval of non-negative consecutive integers.
assert(isvector(N));
assert(all(N >= 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

% Compute dimension.
dim = N(end)^2 + 2*N(end) - N(1)^2 + 1;

% Compute offset.
offset = (N(1)-1)^2 + 2*(N(1)-1);

% Initialise vectors of eigenvalues.
d = zeros(dim, 1);
for k=N
    % Create indices.
    idx = k^2 - offset - 1;
    % Save eigenvalues.
    d(idx+1:idx+2*k+1) = k*(k+1);
end

end