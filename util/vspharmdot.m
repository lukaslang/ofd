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
function Z = vspharmdot(v, F, V, N)
%VSPHARMDOT Computes the dot product between a vector field and vector
%spherical harmonics.
%
%   Z = VSPHARMDOT(v, F, V, N) takes a vector field v on the faces F of a
%   triangulation F, V and computes the dot product with vector spherical
%   harmonics of degrees N.
%
%   Z are scalar fields on the faces and is of size [size(F, 1), dim],
%   where dim is the dimension induced by degrees N.
%
%   Note that parallelisation was chosen this way since for each k the
%   number of dot products to calculate is 2*(2*k+1). Each of the 2*k+1
%   computations can be done with no communcation overhead.

assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(size(v, 2) == 3);
assert(size(v, 1) == size(F, 1));

% Check if N is an interval of positive consecutive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

% Compute dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

% Compute offset.
offset = (N(1)-1)^2 + 2*(N(1)-1);

% Compute dot products.
Z = zeros(size(F, 1), dim);
for k=N
    % Create vector spherical harmonics of degree k.
    [Y1, Y2] = vspharm(k, F, V);
    % Compute index.
    idx = k^2 - offset - 1;
    % Run through all orders.
    parfor l=1:2*k+1
        Z(:, idx + l) = dot(v, squeeze(Y1(:, l, :)), 2);
    end
    % Create indices.
    idx = idx + dim/2;
    parfor l=1:2*k+1
        Z(:, idx + l) = dot(v, squeeze(Y2(:, l, :)), 2);
    end
end
end