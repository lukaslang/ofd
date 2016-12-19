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
function tests = sphericalIntegralTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Create parametrisation and coordinates.
m = 500;
n = 500;

% Create parametrisation and polar coordinates.
[phi, t] = ndgrid(linspace(-pi, pi, m), linspace(-1, 1, n));

for d=0:5
    % Create spherical harmonics.
    N = d;
    Ynj = spharmp(N, phi, t);
    verifyFalse(testCase, isempty(Ynj));
    verifyEqual(testCase, size(Ynj), [m*n, 2*N + 1]);

    for k=1:2*N+1
        fprintf('Checking norm of degree %i, order %i...\n', d, k);
        % Create data.
        f = reshape(Ynj(:, k), m, n).^2;

        % Compute surface integral.
        v = sphericalIntegral(phi, t, f);

        verifyTrue(testCase, isscalar(v));
        verifyEqual(testCase, sqrt(v), 1, 'AbsTol', 1e-3);
    end
end
end

function constantFunctionTest(testCase)

% Create parametrisation and coordinates.
m = 100;
n = 100;

% Create parametrisation and polar coordinates.
[phi, t] = ndgrid(linspace(-pi, pi, m), linspace(-1, 1, n));

% Compute surface integral over function which is constant one.
v = sphericalIntegral(phi, t, ones(m, n));

verifyTrue(testCase, isscalar(v));
verifyEqual(testCase, v, 4*pi, 'AbsTol', 1e-6);

end