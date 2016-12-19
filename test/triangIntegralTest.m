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
function tests = triangIntegralTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function constantFunctionTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
n = size(F, 1);

% Compute surface integral over function which is constant one.
v = triangIntegral(F, V, ones(n, 1));

verifyTrue(testCase, isscalar(v));
verifyEqual(testCase, v, 4*pi, 'AbsTol', 2e-2);

end

function constantFunctionWithAreaTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
n = size(F, 1);

% Compute surface integral over function which is constant one.
v = triangIntegral(F, V, ones(n, 1), triangArea(F, V));

verifyTrue(testCase, isscalar(v));
verifyEqual(testCase, v, 4*pi, 'AbsTol', 2e-2);

end

function vspharmOnbTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);

% Compute triangle areas.
a = triangArea(F, V);

% Check if vector spherical harmonics form an ONB of norm one.
for N=1:5
    % Create vector spherical harmonics.
    [Y1, Y2] = vspharm(N, F, V);

    for k=1:2*N+1
        fprintf('Checking norm of degree %i, order %i...\n', N, k);
        
        % Compute surface integral.
        v = triangIntegral(F, V, dot(Y1(:, k, :), Y1(:, k, :), 3), a);
        verifyTrue(testCase, isscalar(v));
        verifyEqual(testCase, sqrt(v), 1, 'AbsTol', 1e-2);
        
        % Compute surface integral.
        v = triangIntegral(F, V, dot(Y2(:, k, :), Y2(:, k, :), 3), a);
        verifyTrue(testCase, isscalar(v));
        verifyEqual(testCase, sqrt(v), 1, 'AbsTol', 1e-2);
    end
end
end

function vspharmIntegralTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);

% Compute triangle areas.
a = triangArea(F, V);

% Check if surface integral over vector spherical harmonics are zero/one.
for N=1:4
    % Create vector spherical harmonics.
    [Y1, Y2] = vspharm(N, F, V);
    for k=1:2*N+1
        for l=1:2*N+1
            fprintf('Checking int Y_nj^i cdot Y_nk^j of degree n=%i, orders j=%i and k=%i...\n', N, k, l);
            % Compute surface integral.
            v = triangIntegral(F, V, dot(Y1(:, k, :), Y1(:, l, :), 3), a);
            verifyTrue(testCase, isscalar(v));
            if(l == k)
                verifyEqual(testCase, v, 1, 'AbsTol', 1e-2);
            else
                verifyEqual(testCase, v, 0, 'AbsTol', 1e-2);
            end
            % Compute surface integral.
            v = triangIntegral(F, V, dot(Y1(:, k, :), Y2(:, l, :), 3), a);
            verifyTrue(testCase, isscalar(v));
            verifyEqual(testCase, sqrt(v), 0, 'AbsTol', 1e-2);
            % Compute surface integral.
            v = triangIntegral(F, V, dot(Y2(:, k, :), Y2(:, l, :), 3), a);
            verifyTrue(testCase, isscalar(v));
            if(l == k)
                verifyEqual(testCase, v, 1, 'AbsTol', 1e-2);
            else
                verifyEqual(testCase, v, 0, 'AbsTol', 1e-2);
            end
        end
    end
end
end