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
function [a, Ynj] = spharmFit(n, phi, t, f)
%SPHARMFIT Fits spherical harmonics.
%
%   [a, Ynj] = SPHARMFIT(n, phi, t, f) fits spherical harmonics of degree n 
%   through data f at points [phi, t, 1] and returns coefficients a_nj.
%   Note that n is scalar and phi, t, and f are vectors. Moreover, phi is 
%   in [0, 2pi[ and t in [-1, 1].

m = length(phi);

assert(n >= 0);
assert(size(phi, 1) == m);
assert(size(phi, 2) == 1);
assert(size(t, 1) == m);
assert(size(t, 2) == 1);
assert(size(f, 1) == m);
assert(size(f, 2) == 1);

% Create m times (n+1)^2 matrix of spherical harmonics.
Ynj = cell2mat(arrayfun(@(m) spharm(m, phi, t), 0:n, 'UniformOutput', false));

% Solve linear system.
a = cgs(Ynj'*Ynj, Ynj'*f, 10e-6, 30);

end