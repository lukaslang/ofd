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
function U = ofhd(N, F, V, f1, f2, h, alpha, s)
%OFHD Returns a hierarchical optical flow decomposition on a triangulation.
%
%   U = OFHD(N, F, V, f1, f2, h, alpha, s) takes a triangulation F, V and 
%   images f1, f2 on the vertices of the triangulation and returns a 
%   hierarchical optical flow decomposition. The number of components is
%   specified by the length of the vectors alpha and s. Scalar h is a 
%   spacing parameter.
%
%   U is defined on the faces F and is of size [n, 3, k], where
%   n = size(F, 1) is the number of faces and k = length(alpha) = length(s)
%   is the number of components. In the i-th iteration, the regularisation
%   term is defined by alpha(i) and s(i).

m = size(V, 1);
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(isvector(f1));
assert(isvector(f2));
assert(size(f1, 1) == m);
assert(size(f2, 1) == m);
assert(isvector(alpha));
assert(all(alpha >= 0));
assert(isvector(s));
assert(length(alpha) == length(s));
assert(h > 0);
assert(N > 0);

% Identify number of components.
nd = length(alpha);

% Compute functions needed for solving the linear system.
[dim, U, d, b] = linearsystem(F, V, 1:N, f1, f2, h, 1e-6);

% Initialise coefficients matrix.
ud = zeros(dim, nd);

% Decompose into components.
for k=1:nd
    % Solve linear system.
    [u, ~] = ofsolve(dim, U, b, d, alpha(k), s(k), 30);
    % Store coefficients.
    ud(:, k) = u;
    % Adjust right hand side.
    b = b - U*u;
end

% Compute cummulative sum of coefficients.
ud = cumsum(ud, 2);

% Recover functions.
[U1, U2] = vspharmsynth(1:N, F, V, ud);

% Add up vector fields of both types.
U = U1 + U2;

end