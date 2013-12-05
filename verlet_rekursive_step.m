function [ u, v, q, x ] = verlet_rekursive_step( f, u, v, q, h0, n, ...
    depth, max_depth, x, steptype, step_plot )
%VERLET_REKURSIVE_STEP Recursive step for St?rmer-Verlet composition method
%   
%   Parameters:
%       f: Right hand site of the ODE
%       u: Matrix containing results of u
%       v: Matrix containing results of v
%       q: Matrix containing values of u_1/2
%       h0: Original step size (which is splitted into five steps)
%       n: Position in result matrices
%       depth: Current depth
%       x: Absolute step position
%       steptype: Step type (triple or suzuki)
%       step_plot: If set to true plots each step size
%
%   Returns:
%       u: Matrix containing results of u with new solutions added
%       v: Matrix containing results of v with new solutions added
%       q: Matrix containing values of u with new values added
%       x: Absolute step position after evaluation
%
%   Copyright 2013 by Dominik Hagl, Christian Buchmayr

    if strcmp(steptype, 'triple')
        gamma = triple((max_depth-depth)*2+2);
    else
        gamma = suzuki((max_depth-depth)*2+2);
    end

    for t = 1:1:size(gamma,2)
        h = gamma(t) * h0;
        
        if depth == 1
            q(:, n) = u(:, n) + h/2 * f(v(:, n));
            v(:, n+1) = v(:, n) + h*q(:, n);
            u(:, n+1) = q(:, n) + h/2 * f(v(:, n+1));
            
            % visualize step size; see tripple_gamma (note how these 
            % exceed h*[0,1]) and suzuki_gamma (stays in h*[0,1])
            if step_plot
                plot([x, x+h], [n, n+1], '.-b');
            end
        else
            [ u, v, q ] = verlet_rekursive_step( f, u, v, q, h, n, ...
                depth-1, max_depth, x, steptype, step_plot );
        end
        
        x = x + h;
        
        n = n + size(gamma, 2)^(depth-1);
    end
end

