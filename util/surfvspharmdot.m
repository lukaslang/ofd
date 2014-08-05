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
function Z = surfvspharmdot(v, F, V, N, Ns, c)
%VSPHARMDOT Computes the dot product between a vector field on a surface 
%and the pushforward of vector spherical harmonics.
%
%   Z = VSPHARMDOT(v, F, V, N, Ns, c) takes a vector field v on the faces F
%   of a triangulation F, V of the unit spere and computes the dot product 
%   with vector spherical harmonics of degrees N on a surface defined by Ns
%   and c.
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

% Compute synthesis.
[~, rho] = surfsynth(Ns, V, c);

% Compute gradient of rho on triangulation.
g = grad(F, V, rho);

% Compute triangle incenters and project to unit sphere.
TR = TriRep(F, V);
IC = TR.incenters;
len = sqrt(sum(IC .^ 2, 2));
IC = bsxfun(@rdivide, IC, len);
[~, rhoic] = surfsynth(Ns, IC, c);

% Compute dot products.
Z = zeros(size(F, 1), dim);
for k=N
    % Create vector spherical harmonics of degree k.
    [Y1, Y2] = vspharm(k, F, V);
    % Compute pushforward to surface.    
    Y1 = pushforward(Y1, IC, rhoic, g);
    Y2 = pushforward(Y2, IC, rhoic, g);
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