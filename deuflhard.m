function [ q, dq ] = deuflhard(omega, g, q0, dq0, h, N, progress)
%DEUFLHARD ODE solver with Deuflhard's trigonometric method
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
    cosho = cos(ho);
    sincho = sinc(ho/pi);
    sinho = sin(ho);
    
    %initialize second startvalue (two-step formulation)
    q1 = cosho*q0 + sincho*h*dq0 + h^2*sincho*g(q0)/2;
    dq1 = -ho*sinho*q0/h + cosho*dq0 + (h/2)*(g(q1) + cosho*g(q0));
    
    q = [q0 q1 zeros(dim, N-2)];
    dq = [dq0 dq1 zeros(dim, N-2)];

    for n = 3:N;
        q(:, n) = h^2*sincho*g(q(:, n-1)) + ...
            2*cosho*q(:, n-1) - q(:, n-2);
        dq(:, n) = h*(g(q(:, n))-g(q(:, n-2)))/2 + ...
            2*cosho*dq(:, n-1) - dq(:, n-2);
        
        if progress && mod(n, percent) == 0
            fprintf('\b\b\b\b\b%4d%%', n/percent);
        end
    end
    
    if progress
        fprintf('\b\b\b\b\b 100%%\n');
    end
end

