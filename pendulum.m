%% Pendulum from HLW2006, pp 5-7 and example p 188

clear all;

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

% range_p = [1.35 2.5];
% range_q = [-.5 .5];
% h = pi/3;
% N = 3;

range_p = [.35 1.5];
range_q = [5.8 6.8];
h = pi/9;
N = 21;

%[cat, colormap, alpha] = imread('cat24.png');
%image([-.5 .5], [2.5 1.35], cat);

[X, p_range, q_range, colorMap, alpha] = ...
    transform_image('images/cat24.png', range_p, range_q, ...
        @(dp, dq, p0, q0, h, N) verlet(@(q) dp(0, q), p0, q0, h, N), ...
        dp, dq, h, N);

for i = 1:N
    image(q_range{i}, p_range{i}, X{i}, 'AlphaData', alpha{i});
end
%colormap('gray');
colormap(colorMap);

axis([-2 9 -3 3])
