function [u, v] = seuler(du, dv, u0, v0, h, N)
%SEULER Symplectic euler
%   Multidimensional symplectic euler for two functions u and v
%   
%   Parameters
%       du: derivative du/dt
%       dv: derivative dv/dt
%       u0: initial value for u
%       v0: initial value for v
%       h: step size
%       N: number of steps
%
%   Returns
%       u: values of u
%       v: values of v
%
%   Author: Dominik Hagl
    
    fprintf('Symplectic Euler:   0%%');

    dim = length(u0);
    percent = floor( N / 100 );

    u = [u0 zeros(dim, N)];
    v = [v0 zeros(dim, N)];

    for n = 1:N+1;
        % u(n+1) = u(n) + h*du(u(n+1), v(n));
        for i = 1:dim
            e = zeros(1, dim);
            e(i) = 1;
            u(i, n+1) = fzero( ...
                @(g) g - u(i, n) - h * e*feval(du, g, v(:, n)), ...
                u(i, n) ...
            );
        end

        v(:, n+1) = v(:, n) + h*dv(u(:, n+1), v(:, n));

        if mod(n, percent) == 0
            fprintf('\b\b\b\b\b%4d%%', n/percent);
        end
    end
    
    fprintf('\b\b\b\b\b 100%%\n');
end

