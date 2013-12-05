function [] = pendulumH()
%% Pendulum from HLW2006, example p 367

close all;
clear all;

% Problem from p 5
dp = @(p, q) -sin(q);
dq = @(p, q) p;

N = 1000;
h = 1.55;

%Hmod = @(p,q) p.^2/2-cos(q)+h^2/48*(cos(2*q)-2*p.^2.*cos(q));
Hamiltonian = @(p,q) p.^2/2-cos(q);

LineStyle = {'k-', 'r-', 'b-'};

figure;
hold on;
grid on;

%% Solve ODE with multiple initial values

Inits = [0, -1.5; 0, -2.5; 1.5, pi];
b = size(Inits);
for i = 1:b(1)
    % St?rmer-Verlet scheme
    [p, q] = verlet(@(q) dp(0, q), Inits(i,1), Inits(i,2), h, N);

    % Plot
    H = Hamiltonian(p,q);

    p = plot(0:N, H, char(LineStyle(i)));
    set(p,'LineWidth',1.5);
end

legend(getLegend(Inits(:,1), Inits(:,2), '(p_0,q_0)=(%1.2f, %1.2f)'));

% ylim([0.5,3]);
end

function s = getLegend(values, values2, format)
    s = {};

    for k = 1:length(values)
        s(k) = cellstr(sprintf(format, values(k), values2(k)));
    end
end

