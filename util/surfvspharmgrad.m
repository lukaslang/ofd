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
function [Z1, Z2, Z3, Z4] = surfvspharmgrad(X1, X2, F, V, N, Ns, c)
%SURFVSPHARMGRAD Computes covariant derivative of vector spherical harmonics.
%
%   [Z1, Z2, Z3, Z4] = SURFVSPHARMGRAD(X1, X2, F, V, N, Ns, c) takes vector
%   fields X1, X2 defined on the faces F of a triangulation F, V of the unit
%   sphere such that X1, X2 form a basis for the tangent space. A sphere-like
%   surface specified by Ns and c, and computes the covariant derivative of
%   the pushforward of vector spherical harmonics of degrees N with respect
%   to the vector fields X1 and X2.
%
%   Note that X1 and X2 are of size [size(F, 1), 3].
%
%   Zi are scalars on the faces and are of size [size(F, 1), dim], where 
%   dim is the dimension induced by degrees N.
%
%   Note that parallelisation was chosen this way since for each degree k
%   the number of dot products to calculate is 2*(2*k+1). Each of the 2*k+1
%   computations can be done with no communcation overhead.

assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(size(X1, 1) == size(F, 1));
assert(size(X1, 2) == 3);
assert(size(X2, 1) == size(F, 1));
assert(size(X2, 2) == 3);

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
[Vs, rho] = surfsynth(Ns, V, c);

% Compute gradient of rho on triangulation.
g = grad(F, V, rho);

% Compute triangulation.
T = TriRep(F, Vs);
DV = T.incenters;
% Compute triangles spanned by incenters of neighboring faces.
DF = T.neighbors;
% Compute triangle heights.
H = height(DF, DV);
% Compute squared lengths.
lenH = sum(H(:, 2:3, :).^2, 3);

% Compute face normals of triangulation of incenters.
DT = TriRep(DF, DV);
FN = -DT.faceNormals;

% Compute triangle incenters and project to unit sphere.
TR = TriRep(F, V);
IC = TR.incenters;
len = sqrt(sum(IC .^ 2, 2));
IC = bsxfun(@rdivide, IC, len);
[~, rhoic] = surfsynth(Ns, IC, c);

% Compute covariant derivatives.
Z1 = zeros(size(F, 1), dim);
Z2 = zeros(size(F, 1), dim);
Z3 = zeros(size(F, 1), dim);
Z4 = zeros(size(F, 1), dim);
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
        % Compute covariant derivatives of Y1.
        D1 = vgrad(DF, DV, squeeze(Y1(:, l, :)), X1, H, lenH, FN);
        D2 = vgrad(DF, DV, squeeze(Y1(:, l, :)), X2, H, lenH, FN);
        % Compute components.
        Z1(:, idx + l) = dot(D1, X1, 2);
        Z2(:, idx + l) = dot(D1, X2, 2);
        Z3(:, idx + l) = dot(D2, X1, 2);
        Z4(:, idx + l) = dot(D2, X2, 2);
    end
    % Create indices.
    idx = idx + dim/2;
    parfor l=1:2*k+1
        % Compute covariant derivatives of Y2.
        D1 = vgrad(DF, DV, squeeze(Y2(:, l, :)), X1, H, lenH, FN);
        D2 = vgrad(DF, DV, squeeze(Y2(:, l, :)), X2, H, lenH, FN);
        % Compute components.
        Z1(:, idx + l) = dot(D1, X1, 2);
        Z2(:, idx + l) = dot(D1, X2, 2);
        Z3(:, idx + l) = dot(D2, X1, 2);
        Z4(:, idx + l) = dot(D2, X2, 2);
    end
end
end