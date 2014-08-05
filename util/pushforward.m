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
function Z = pushforward(Y, X, rho, g)
%PUSHFORWARD Computes the pushforward of tangent vector fields on the 
%unit sphere to a sphere-like surface.
%
%   Z = PUSHFORWARD(Y, X, rho, g) takes vector fields Y defined at points X
%   on the unit sphere, a scalar rho defined at X and the gradient g at X,
%   and returns the pushforward Z of Y to the sphere-like surface defined 
%   by rho and g.
%
%   Note that X is of size m-by-3 and Y is an m-by-k-by-3 matrix defining 
%   vectors at X on the unit sphere. rho is a vector of length m and g is a
%   matrix of size m-by-3.
%
%   Z is a matrix of size m-by-k-by-3.

% Check if dimensions comply.
assert(isvector(rho));
assert(size(rho, 1) == size(Y, 1));
assert(size(Y, 1) == size(X, 1));
assert(size(Y, 3) == 3);
assert(size(X, 2) == 3);
assert(size(g, 1) == size(Y, 1));
assert(size(g, 2) == 3);

% Compute pushforward for each vector field.
Z = zeros(size(Y));
parfor k=1:size(Y, 2)
    Yk = squeeze(Y(:, k, :));
    Z(:, k, :) = bsxfun(@times, Yk, rho) + bsxfun(@times, X, sum(g .* Yk, 2));
end