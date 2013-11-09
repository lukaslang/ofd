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
function [a, Ynj] = spharmpFit(N, phi, t, f)
%SPHARMPFIT Fits spherical harmonics using a cylindrical parametrisation.
%
%   [a, Ynj] = SPHARMPFIT(n, phi, t, f) fits spherical harmonics of degree N 
%   through data f at grid points phi, t and returns coefficients a_nj.
%   Note that N is scalar and phi, t, and f are m-times-n matrices. 
%   Moreover, phi is in [0, 2pi[ and t in [-1, 1].

[m, n] = size(phi);

assert(N >= 0);
assert(size(phi, 1) == m);
assert(size(phi, 2) == n);
assert(size(t, 1) == m);
assert(size(t, 2) == n);
assert(size(f, 1) == m);
assert(size(f, 2) == n);

% Create m times (n+1)^2 matrix of spherical harmonics.
Ynj = cell2mat(arrayfun(@(m) spharmp(m, phi, t), 0:N, 'UniformOutput', false));

% Solve linear system.
a = cgs(Ynj'*Ynj, Ynj'*f(:), 10e-6, 30);

end