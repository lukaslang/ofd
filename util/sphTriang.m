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
function [F, V] = sphTriang(n)
%SPHTRIANG Computes a triangulation of the unit sphere.
%
%   [F, V] = SPHTRIANG returns faces F and vertices V of a triangulation of
%   the unit sphere generated from the icosahedron. F is the 12x3 matrix 
%   with the row indices of the 3 points in the vertex matrix V. Note that 
%   indices are in clockwise order. V itself is a 20x3 matrix with points 
%   in R3.
%
%   [F, V] = SPHTRIANG(n) does n > 0 refinements of the mesh.
%
%   Note that vertices lie on the unit sphere.
%
%   Usage with TriRep:
%   [F, V] = sphTriang;
%   DT = TriRep(F, V);

if(nargin == 0)
    n = 0;
end

% Create icosahedron.
[F, V] = icosahedron;

% Project vertices to unit sphere.
V = normalise(V);

% Refine mesh.
for k=1:n
    numV = size(V, 1);
    % Refine mesh.
    [F, V] = refineMesh(F, V);
    % Normalise new vertices.
    V(numV+1:end, :) = normalise(V(numV+1:end, :));
end

end

function V = normalise(V)
%NORMALISE Normalises row vectors.
%
%   V = normalise(V) takes a n-by-3 matrix V and returns normalises row
%   vectors.

len = sqrt(sum(V.^2, 2));
V = bsxfun(@rdivide, V, len);

end

function [F, V] = icosahedron

% Compute golden ratio.
phi = (1 + sqrt(5))/2;

% Create icosahedron from cartesian coordinates.
V = [   0,    1,  phi;
        0,    1, -phi;
        0,   -1,  phi;
        0,   -1, -phi;
        1,  phi,    0;
        1, -phi,    0;
       -1,  phi,    0;
       -1, -phi,    0;
      phi,    0,    1;
      phi,    0,   -1;
     -phi,    0,    1;
     -phi,    0,   -1];
 
% Create faces in clockwise order.
F = [   1,  3, 11;
        1, 11,  7;
        1,  7,  5;
        1,  5,  9;
        1,  9,  3;
        3,  6,  8;
        3,  9,  6;
        3,  8, 11;
        9,  5, 10;
        9, 10,  6;
       11,  8, 12;
       11, 12,  7;
        8,  6,  4;
        6, 10,  4;
        8,  4, 12;
        5,  2, 10;
        5,  7,  2;
        7, 12,  2;
        2, 12,  4;
        4, 10,  2];
end