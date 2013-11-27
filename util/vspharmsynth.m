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
function U = vspharmsynth(N, F, V, u)
%VSPHARMSYNTH Computes vector spherical harmonics synthesis.
%
%   U = VSPHARMSYNTH(N, F, V, u) takes coefficients u for vector sphercial
%   harmonics of degrees N and a triangulation F, V and returns a vector 
%   field U defined on the faces F.
%
%   U is of size [size(F, 1), 3].
%
%   Note that N is a vector of positive consecutive integers!

% Compute and check dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);
assert(isvector(u));
assert(length(u) == dim);

% Check if N is an interval of consecutive positive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

% Get number of faces.
n = size(F, 1);

% Compute offset for interval.
offset = (N(1)-1)^2 + 2*(N(1)-1);

% Recover vector field.
U = zeros(n, 3);
parfor k=N
    % Create vector spherical harmonics of degree k.
    [Y1, Y2] = vspharm(k, F, V);
    % Create temporary variables with coefficients.
    idx = [k^2, k^2+2*k] - offset;
    u1 = u(idx(1):idx(2));
    u2 = u(idx(1)+dim/2:idx(2)+dim/2);
    % Compute vectors.
    U = U + cat(2, sum(bsxfun(@times, squeeze(Y1(:, :, 1)), u1'), 2), sum(bsxfun(@times, squeeze(Y1(:, :, 2)), u1'), 2), sum(bsxfun(@times, squeeze(Y1(:, :, 3)), u1'), 2));
    U = U + cat(2, sum(bsxfun(@times, squeeze(Y2(:, :, 1)), u2'), 2), sum(bsxfun(@times, squeeze(Y2(:, :, 2)), u2'), 2), sum(bsxfun(@times, squeeze(Y2(:, :, 3)), u2'), 2));
end
end