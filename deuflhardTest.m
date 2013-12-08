
%% Problem from p 17
omega = 50;
A = -omega^2/2.*[1 -1;-1 1];
g = @(q) [4*(q(1)-q(2))^3; 4*(q(1)-q(2))^3] ;
ddq = @(q) A*q + g(q); 

% initial values
q0 = [0;1];
dq0 = [1;0];


%% Step size and number of steps
N = 100;    % 1e4 steps
h = 0.001; 

%% Solve ODE
tic

% Deuflhard method
[q, dq] = deuflhard(A^0.5, g, q0, dq0, h, N, true);

toc
