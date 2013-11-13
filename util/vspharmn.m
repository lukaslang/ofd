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
function [Y, d] = vspharmn(N, F, V)
%VSPHARMN Creates vector spherical harmonics up to a certain degree.
%
%   [Y, D] = VSPHARMN(N, F, V) takes a degree N > 0 and a triangulation F,
%   V and computes vector spherical harmonics Y_nj^i of degree n=1,...N and
%   order j=-n,...,n for every face F. Matrix Y = [Y1, Y2] is composed of
%   vector spherical harmonics of type 1 and 2. The vector d contains the
%   eigenvalues n*(n+1).
%
%   Note that the matrix Y contains vectors in R3 and is of size 
%   f-by-g-by-3, where f = size(F, 1) is the number of faces, 
%   g = 2*(N^2 + 2*N) is the dimension of the set of vector spherical 
%   harmonics up to degree N.
%
%   Vector d is of length 2*(N^2 + 2*N).

n = size(F, 1);

% Create vector containing eigenvalues.
d = zeros(N^2 + 2*N, 1);

% Create vector spherical harmonics.
Y1 = zeros(n, N^2 + 2*N, 3);
Y2 = zeros(n, N^2 + 2*N, 3);

c = 1;
for k=1:N
    % Generate vector spherical harmonics of degree k.
    [Yi, Yj] = vspharm(k, F, V);
    Y1(:, c:(c+2*k), :) = Yi;
    Y2(:, c:(c+2*k), :) = Yj;
    % Save eigenvalues.
    d(c:(c+2*k)) = repmat(k*(k+1), (c+2*k)-c+1, 1);
    c = c + 2*k + 1;
end
Y = cat(2, Y1, Y2);
d = repmat(d, 2, 1);

end