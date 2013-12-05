%% Pendulum from HLW2006, pp 5-7 and example p 188

% Problem from p 5
dp = @(p, q) -sin(q);
dq = @(p, q) p;

N = 100;
h = .1;

figure;
hold on;
grid on;

%% Solve ODE with multiple initial values

for q0 = -2*pi:pi/2:4*pi
    for p0 = -3:.5:3
        % St?rmer-Verlet scheme
        [p, q] = verlet(@(q) dp(0, q), p0, q0, h, N);

        % Plot
        plot(q, p);
    end
end

%% Plot area preserved shape
% Q0 = -1:.1:0;
% P0 = 1.5:.1:2.5;
% h = pi/3;
% N = 3;

Q0 = 5.8:.1:6.8;
P0 = .5:.1:1.5;
h = pi/3;
N = 8;

n = length(Q0) * length(P0);
q = cell(1, N);
p = cell(1, N);

for i = 1:N
    q{i} = zeros(1, n);
    p{i} = zeros(1, n);
end

j = 0;

for q0 = Q0
    for p0 = P0
        j = j+1;
        
        [p_temp, q_temp] = verlet(@(q) dp(0, q), p0, q0, h, N);
        
        for i = 1:N
            q{i}(j) = q_temp(i);
            p{i}(j) = p_temp(i);
        end
    end
end

%% FIXME: Just plot border and connect correct neighbors

colorOrder = get(gca,'ColorOrder');
[m, n] = size(colorOrder);
for i = 1:N
    plot(q{i}, p{i}, '.-', 'LineWidth', 5, 'MarkerSize', 5, ...
        'Color', colorOrder(mod(i, m)+1, :));
end

axis([-2 9 -3 3])
