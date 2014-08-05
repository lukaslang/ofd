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
function U = ofs(N, Ns, c, F, V, f1, f2, h, alpha)
%OFS Computes the optical flow on a sphere-like surface.
%
%   U = OFS(N, Ns, c, F, V, f1, f2, h, alpha) takes a triangulation F, V of
%   a surface specified by Ns and C, and images f1, f2 on the vertices of 
%   the triangulation and returns the optical flow U. Scalar h is a spacing
%   parameter and alpha is the regularisation parameter.
%
%   U is defined on the faces F and is of size [size(F, 1), 3].

m = size(V, 1);
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(isvector(f1));
assert(isvector(f2));
assert(size(f1, 1) == m);
assert(size(f2, 1) == m);
assert(alpha >= 0);
assert(h > 0);
assert(N > 0);
assert(isscalar(N));

% TODO: add assertions for Ns and c!
assert(isvector(c));

% Compute linear system for optical flow.
[~, A, D, b] = surflinearsystem(F, V, Ns, c, 1:N, f1, f2, h, 1e-6);

% Solve linear system.
u = gmres(A + alpha * D, b, [], 1e-6, 30);

% Recover function.
[U1, U2] = vspharmsynth(1:N, F, V, u);
U(:, 1, :) = U1 + U2;

% Pushforward vector field to surface.
U = squeeze(surfpushforward(Ns, c, F, V, U));

end