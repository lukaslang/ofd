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
function [Y1, Y2] = vspharmp(N, phi, t)
%VSPHARMP Generates fully normalised vector spherical harmonics using a
%cylindrical parametrisation.
%
%   [Y1, Y2] = VSPHARMP(N, phi, t) takes grid matrices phi, t of size 
%   m-times-n (e.g. obtained by ndgrid) and returns vector spherical 
%   harmonics Y_Nj of degree N > 0 and j=-N,...,N. phi must be in [0, 2pi[
%   and t in [-1, 1].
%
%   Note that size(Yi) = [m, n, 2*N + 1, 3] for i={1, 2}.

[m, n] = size(phi);

assert(N > 0);
assert(size(phi, 1) == m);
assert(size(phi, 2) == n);
assert(size(t, 1) == m);
assert(size(t, 2) == n);

% Compute spherical harmonics.
Ynj = spharmp(N, phi, t);

[dphi, dt] = gradient(reshape(Ynj, m, n, 2*N + 1));

% Compute orthonormal basis.
ephi(:, :, 1) = -sin(phi);
ephi(:, :, 2) = cos(phi);
ephi(:, :, 3) = 0;

et(:, :, 1) = -t .* cos(phi);
et(:, :, 2) = -t .* sin(phi);
et(:, :, 3) = sqrt(1 - t.^2);

er(:, :, 1) = sqrt(1 - t.^2) .* cos(phi);
er(:, :, 2) = sqrt(1 - t.^2) .* sin(phi);
er(:, :, 3) = t;

% Compute vectorial spherical harmonics.
Y1(:, :, :, 1) = repmat(ephi(:, :, 1) ./ sqrt(1 - t.^2), [1, 1, 2*N + 1]) .* dphi + repmat(sqrt(1 - t.^2) .* et(:, :, 1), [1, 1, 2*N + 1]) .* dt;
Y1(:, :, :, 2) = repmat(ephi(:, :, 2) ./ sqrt(1 - t.^2), [1, 1, 2*N + 1]) .* dphi + repmat(sqrt(1 - t.^2) .* et(:, :, 2), [1, 1, 2*N + 1]) .* dt;
Y1(:, :, :, 3) = repmat(ephi(:, :, 3) ./ sqrt(1 - t.^2), [1, 1, 2*N + 1]) .* dphi + repmat(sqrt(1 - t.^2) .* et(:, :, 3), [1, 1, 2*N + 1]) .* dt;

Y2(:, :, :, 1) = Y1(:, :, :, 2) .* repmat(er(:, :, 3), [1, 1, 2*N + 1]) - Y1(:, :, :, 3) .* repmat(er(:, :, 2), [1, 1, 2*N + 1]);
Y2(:, :, :, 2) = Y1(:, :, :, 3) .* repmat(er(:, :, 1), [1, 1, 2*N + 1]) - Y1(:, :, :, 1) .* repmat(er(:, :, 3), [1, 1, 2*N + 1]);
Y2(:, :, :, 3) = Y1(:, :, :, 1) .* repmat(er(:, :, 2), [1, 1, 2*N + 1]) - Y1(:, :, :, 2) .* repmat(er(:, :, 1), [1, 1, 2*N + 1]);

% Normalise.
Y1 = Y1 ./ sqrt(N*(N+1));
Y2 = Y2 ./ sqrt(N*(N+1));

end