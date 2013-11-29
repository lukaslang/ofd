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
function U = vspharmsynth(N, F, V, u, mem)
%VSPHARMSYNTH Computes vector spherical harmonics synthesis.
%
%   U = VSPHARMSYNTH(N, F, V, u) takes coefficients u for vector sphercial
%   harmonics of degrees N and a triangulation F, V and returns a vector 
%   field U defined on the faces F. 
%
%   Note that u can be a matrix of size [size(F, 1), n] where columns are 
%   the coefficients.
%
%   U = VSPHARMSYNTH(N, F, V, u, mem) additionally takes a memory
%   constraint in bytes and allows VSPHARMSYNTH to use up to mem bytes.
%
%   Note that if mem is specified, VSPHARMSYNTH then creates several 
%   degrees of vector spherical harmonics at once using mem bytes. However,
%   at least one degree is generated possibly exceeding the specified
%   memory!
%
%   If u is a vector, then U is of size [size(F, 1), 3]. If u is a matrix 
%   of size [size(F, 1), m] then U is of size [size(F, 1), 3, m].
%
%   Note that N is a vector of positive consecutive integers!

% Compute and check dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);
assert(size(u, 1) == dim);

% Check if N is an interval of consecutive positive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

% Get number of coefficient vectors.
m = size(u, 2);

% Get number of faces.
n = size(F, 1);

% Compute offset for interval.
offset = (N(1)-1)^2 + 2*(N(1)-1);

% Compute intervals I that fit into mem.
I = {};
if(nargin == 5)
    assert(mem > 0);
    s = N(1);
    while(s <= N(end))
        e = interval(N, s, n, mem);
        I{end+1} = s:e;
        s = e + 1;
    end
else
    I = mat2cell(N', ones(length(N), 1));
end
    
% Create several degrees at once.
U = zeros(n, 3, m);
% Note that for the memory constraints this loop must not run in parallel
% since it would otherwise use more resources! However, the inner loop
% might be subject to parallelisation.
for k=1:length(I)
    % Generate interval.
    M = I{k};
    % Create vector spherical harmonics for interval M.
    [Y1, Y2] = vspharmn(M, F, V);
    % Create indices.
    idx = [M(1)^2, M(end)^2+2*M(end)] - offset;
    % Compute for each column in u.
    for c=1:m
        % Create temporary variables with coefficients.
        u1 = u(idx(1):idx(2), c);
        u2 = u(idx(1)+dim/2:idx(2)+dim/2, c);
        % Compute vectors.
        U(:, :, c) = U(:, :, c) + cat(2, sum(bsxfun(@times, squeeze(Y1(:, :, 1)), u1'), 2), sum(bsxfun(@times, squeeze(Y1(:, :, 2)), u1'), 2), sum(bsxfun(@times, squeeze(Y1(:, :, 3)), u1'), 2));
        U(:, :, c) = U(:, :, c) + cat(2, sum(bsxfun(@times, squeeze(Y2(:, :, 1)), u2'), 2), sum(bsxfun(@times, squeeze(Y2(:, :, 2)), u2'), 2), sum(bsxfun(@times, squeeze(Y2(:, :, 3)), u2'), 2));
    end
    % Explicitly free memory.
    Y1 = [];
    Y2 = [];
end

end

function e = interval(N, s, n, mem)
    % Compute interval ending such that
    % 2*(e^2 + 2*e^2 - s^2 + 1)*n*3*8 = mem.
    e = min(N(end), floor(sqrt(s^2 + mem/(2*n*3*8)) - 1));
end