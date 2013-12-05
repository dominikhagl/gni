function [ u, v ] = verlet(f, u0, v0, h, N, progress)
%VERLET ODE solver with the St?rmer-Verlet scheme
%   Evaluates a ODE of the form
%       du/dt = v
%       dv/dt = f(v)
%   with starting values u0 for u and v0 for v at t=0
%   Function runs N cycles with step size h.
%   
%   Parameters:
%       f: Right hand site of the ODE
%       u0: Start value of u
%       v0: Start value of v
%       h: Step size
%       N: Maximum number of steps
%
%   Returns:
%       u: Approximated values of u at N points
%       v: Approximated values of v at N points
%   
%   See HLW2006, pp 8
%
%   Copyright 2013 by Dominik Hagl, Christian Buchmayr

    if nargin < 6
        progress = false;
    end

    dim = length(u0);
    
    if progress
        fprintf('St?rmer-Verlet:   0%%');
        percent = floor( N / 100 );
    end

    u = [u0 zeros(dim, N)];
    v = [v0 zeros(dim, N)];
    q = zeros(dim, N+1); % u_1/2

    for n = 1:N;
        q(:, n) = u(:, n) + h/2 * f(v(:, n));
        v(:, n+1) = v(:, n) + h*q(:, n);
        u(:, n+1) = q(:, n) + h/2 * f(v(:, n+1));

        if progress && mod(n, percent) == 0
            fprintf('\b\b\b\b\b%4d%%', n/percent);
        end
    end
    
    if progress
        fprintf('\b\b\b\b\b 100%%\n');
    end
end

