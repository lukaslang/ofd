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
function test_suite = triangIntegralTest
    initTestSuite;
end

function constantFunctionTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
n = size(F, 1);

% Compute surface integral over function which is constant one.
v = triangIntegral(F, V, ones(n, 1));

assertTrue(isscalar(v));
assertAlmostEqual(v, 4*pi, 1e-2);

end

function constantFunctionWithAreaTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
n = size(F, 1);

% Compute surface integral over function which is constant one.
v = triangIntegral(F, V, ones(n, 1), triangArea(F, V));

assertTrue(isscalar(v));
assertAlmostEqual(v, 4*pi, 1e-2);

end

function vspharmOnbTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);

% Check if vector spherical harmonics form an ONB of norm one.
for N=1:5
    % Create vector spherical harmonics.
    [Y1, Y2] = vspharm(N, F, V);

    for k=1:2*N+1
        fprintf('Checking norm of degree %i, order %i...\n', N, k);
        
        % Compute surface integral.
        v = triangIntegral(F, V, dot(Y1(:, k, :), Y1(:, k, :), 3));
        assertTrue(isscalar(v));
        assertAlmostEqual(sqrt(v), 1, 1e-2);
        
        % Compute surface integral.
        v = triangIntegral(F, V, dot(Y2(:, k, :), Y2(:, k, :), 3));
        assertTrue(isscalar(v));
        assertAlmostEqual(sqrt(v), 1, 1e-2);
    end
end
end

function vspharmIntegralTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);

% Check if surface integral over vector spherical harmonics are zero/one.
for N=1:4
    % Create vector spherical harmonics.
    [Y1, Y2] = vspharm(N, F, V);
    for k=1:2*N+1
        for l=1:2*N+1
            fprintf('Checking int Y_nj^i cdot Y_nk^j of degree n=%i, orders j=%i and k=%i...\n', N, k, l);
            % Compute surface integral.
            v = triangIntegral(F, V, dot(Y1(:, k, :), Y1(:, l, :), 3));
            assertTrue(isscalar(v));
            if(l == k)
                assertAlmostEqual(v, 1, 1e-2);
            else
                assertAlmostEqual(v, 0, 1e-2);
            end
            % Compute surface integral.
            v = triangIntegral(F, V, dot(Y1(:, k, :), Y2(:, l, :), 3));
            assertTrue(isscalar(v));
            assertAlmostEqual(sqrt(v), 0, 1e-2);
            % Compute surface integral.
            v = triangIntegral(F, V, dot(Y2(:, k, :), Y2(:, l, :), 3));
            assertTrue(isscalar(v));
            if(l == k)
                assertAlmostEqual(v, 1, 1e-2);
            else
                assertAlmostEqual(v, 0, 1e-2);
            end
        end
    end
end
end