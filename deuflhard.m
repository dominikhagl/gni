function [ q, dq ] = deuflhard(omega, g, q0, dq0, h, N, progress)
%VERLET ODE solver with Deuflhard's trigonometric method
%   Evaluates a ODE of the form
%       d^2(q)/dt^2 = -omega^2*q+g(q)
%   with starting values q0 for q and dq0 for d(q)/dt at t=0
%   Function runs N cycles with step size h.
%   
%   Parameters:
%       omega: paramter for quadratic potential
%       g: nonsymmetric right hand side of the ODE
%       q0: Start value of q
%       dq0: Start value of dq
%       h: Step size
%       N: Maximum number of steps
%
%   Returns:
%       q: Approximated values of q at N points
%       dq: Approximated values of d(q)/dt at N points
%   
%   See HLW2006, pp 474
%
%   Copyright 2013 by Dominik Hagl, Christian Buchmayr

    if nargin < 6
        progress = false;
    end

    dim = length(q0);
    
    if progress
        fprintf('Deuflhards Trigonometric Method:   0%%');
        percent = floor( N / 100 );
    end
    
    ho = h*omega;
    
    %initialize second startvalue (two-step formulation)
    q1(:) = cos(ho)*q0(:) + sinc(ho)*h*dq0(:) + h^2*sinc(ho)*g(q0(:));
    dq1(:) = -ho*sin(ho)*q0(:) + cos(ho)*h*dq0(:) + h^2*(g(q1(:)) + cos(ho)*g(q0(:)));
    
    q = [q0 q1 zeros(dim, N)];
    dq = [dq0 dq1 zeros(dim, N)];

    for n = 1:N;
        q(:, n) = h^2*sinc(ho)*g(q(:, n-1)) + 2*cos(ho)*q(:, n-1) - q(:, n-2);
        dq(:, n) = h*(g(:, n)-g(:, n-2))/2 + 2*cos(ho)*dq(:, n-1) - dq(:, n-2);

        if progress && mod(n, percent) == 0
            fprintf('\b\b\b\b\b%4d%%', n/percent);
        end
    end
    
    if progress
        fprintf('\b\b\b\b\b 100%%\n');
    end
end

