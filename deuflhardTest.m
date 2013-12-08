
%% Problem from p 21f
omega = 50;
% A = -omega^2/2.*[1 -1;-1 1];
% g = @(q) [4*(q(1)-q(2))^3; 4*(q(1)-q(2))^3] ;
% ddq = @(q) A*q + g(q); 

A = zeros(6,6);
A(4:6, 4:6) = omega^2*eye(3,3);
g = @(x) -[(x(1)-x(4)).^3-(x(2)-x(5)-x(1)-x(4)).^3;...
            (x(2)-x(5)-x(1)-x(4)).^3-(x(3)-x(6)-x(2)-x(5)).^3;...
            (x(3)-x(6)-x(2)-x(5)).^3+(x(3)+x(6)).^3;...
            -(x(1)-x(4)).^3-(x(2)-x(5)-x(1)-x(4)).^3;...
            -(x(2)-x(5)-x(1)-x(4)).^3-(x(3)-x(6)-x(2)-x(5)).^3;...
            -(x(3)-x(6)-x(2)-x(5)).^3+(x(3)+x(6)).^3];

% initial values
q0 = [1;0;0;omega^(-1);0;0];
dq0 = [1;0;0;1;0;0];


%% Step size and number of steps
N = 25000;    % 1e4 steps
h = .01; 

%% Solve ODE
tic

% Deuflhard method
[x, y] = deuflhard(A^0.5, g, q0, dq0, h, N, true);
% [y, x] = verlet(@(x) -A*x+g(x), dq0, q0, h, N, true);

toc

%% Plot result
I1 = 0.5*(y(4, :).^2 + omega^2*x(4, :).^2);
I2 = 0.5*(y(5, :).^2 + omega^2*x(5, :).^2);
I3 = 0.5*(y(6, :).^2 + omega^2*x(6, :).^2);
I = I1+I2+I3;

figure
plot(1:size(I, 2), [I; I1; I2; I3])

H = @(y, x) .5 * sum(y.^2 + x.^2) + ...
    omega^2 * sum(x(4:6, :).^2) + ...
    .25 * (x(1, :) - x(4, :)).^4 + ...
    .25 * sum(x(2:3, :) - x(5:6, :) - x(1:2, :) - x(4:5, :)).^4 + ...
    .25 * (x(3, :) + x(6, :)).^4;

hold on
plot(1:size(I, 2), H(y, x)-.8);