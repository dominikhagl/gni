
%% Problem from p 17
omega = 50;
% A = -omega^2/2.*[1 -1;-1 1];
% g = @(q) [4*(q(1)-q(2))^3; 4*(q(1)-q(2))^3] ;
% ddq = @(q) A*q + g(q); 

A = zeros(6,6);
A(4:6, 4:6) = -omega^2*eye(3,3);
g = @(x) [(x(1)-x(4)).^3-(x(2)-x(5)-x(1)-x(4)).^3;...
            (x(2)-x(5)-x(1)-x(4)).^3-(x(3)-x(6)-x(2)-x(5)).^3;...
            (x(3)-x(6)-x(2)-x(5)).^3+(x(3)+x(6)).^3;...
            -(x(1)-x(4)).^3-(x(2)-x(5)-x(1)-x(4)).^3;...
            -(x(2)-x(5)-x(1)-x(4)).^3-(x(3)-x(6)-x(2)-x(5)).^3;...
            -(x(3)-x(6)-x(2)-x(5)).^3+(x(3)+x(6)).^3];

% initial values
q0 = [1;0;0;omega^(-1);0;0];
dq0 = [1;0;0;1;0;0];


%% Step size and number of steps
N = 100;    % 1e4 steps
h = 0.001; 

%% Solve ODE
tic

% Deuflhard method
[q, dq] = deuflhard(A^0.5, g, q0, dq0, h, N, true);

toc
