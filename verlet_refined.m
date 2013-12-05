function [ u, v ] = verlet_refined(f, u0, v0, h0, N, max_depth, ...
    steptype, progress, step_plot)
%VERLET_REFINED Composition method of the St?rmer-Verlet scheme
%   Evaluates a ODE of the form
%       du/dt = v
%       dv/dt = f(u)
%   with starting values u0 for u and v0 for v at t=0
%   The method runs N cycles. Each cylce is split into mutiple sub steps
%   calculated by 'triple' or 'suzuki'.
%   
%   Parameters:
%       f: Right hand site of ODE
%       u0: Start value of u
%       v0: Start value of v
%       h0: Initial step size
%       N: Maximum number of steps
%       max_depth: Each step is devided max_depth times
%       steptype: type of steps are used to gain higher order ('triple'
%           or 'suzuki')
%       progress: If set to true progress is displayed
%       step_plot: If set to true the step size is plotted
%
%   Returns:
%       u: Approximated values of u
%       v: Approximated values of v
%   
%   See HLW2006, pp 8, 44-46
%
%   Copyright 2013 by Dominik Hagl, Christian Buchmayr

    if nargin < 8
        progress = false;
    end
    
    if nargin < 9
        step_plot = false;
    end
    
    dim = length(u0);
    percent = floor( N / 100 );
    
    if strcmp(steptype, 'triple')
        step = 3^max_depth;
    else
        steptype = 'suzuki';
        step = 5^max_depth;
    end
    
    if progress
        fprintf('St?rmer-Verlet with %6s jump:   0%%', steptype);
    end
    
    steps = step * N;
    u = [u0 zeros(dim, steps)];
    v = [v0 zeros(dim, steps)];
    q = zeros(dim, steps+1); % u_1/2

    if step_plot
        figure, hold on, grid on;
    end
    
    for n = 1:N;
        [u, v, ~, ~] = verlet_rekursive_step(...
            f, u, v, q, h0, step*(n-1)+1, max_depth, max_depth, ...
            (n-1)*h0, steptype, step_plot);
    
        if mod(n, percent) == 0
            fprintf('\b\b\b\b\b%4d%%', n/percent);
        end
    end
    
    if step_plot
        hold off;
    end
    
    u = u(:, 1:step:end);
    v = v(:, 1:step:end);
    
    fprintf('\b\b\b\b\b 100%%\n');
end

