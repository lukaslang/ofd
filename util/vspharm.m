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
function [Y1, Y2] = vspharm(N, F, V)
%VSPHARM Generates fully normalised vector spherical harmonics on a
%triangulation.
%
%   [Y1, Y2] = VSPHARM(N, F, V) takes a triangulation F, V of the unit 
%   sphere and returns fully normalised vector spherical harmonics Y_Nj of 
%   degree N > 0 and j=-N,...,N.
%
%   Note that size(Yi) = [n, 2*N + 1, 3] for i={1, 2}, where n is the
%   number of faces F.

assert(N > 0);
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);

n = size(F, 1);

% Compute face normals.
T = TriRep(F, V);
Fn = T.faceNormals;

% Compute triangle heigts.
H = height(F, V);

% Compute spherical harmonics.
Ynj = spharm(N, V);

% Compute vector spherical harmonics.
Y1 = zeros(n, 2*N + 1, 3);
Y2 = zeros(n, 2*N + 1, 3);
for k=1:2*N+1
    Y1(:, k, :) = grad(F, V, Ynj(:, k), H);
    Y2(:, k, :) = cross(squeeze(Y1(:, k, :)), Fn);
end

% Normalise.
Y1 = Y1 ./ sqrt(N*(N+1));
Y2 = Y2 ./ sqrt(N*(N+1));

end