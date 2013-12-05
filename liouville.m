%%
% H: Hamiltionian
% 
% M: Number of initial values
% N: Number of steps for ODE evaluation
%
% m: nummer of subplots (rows)
% n: number of subplots (columns)

clear all;

H = @(p, q) p.^2/2 - cos(q);

dp = @(p, q) -sin(q);
dq = @(p, q) p;

M = 15;
N = 200;
h = .1;
delta = 1e-3;

% number of subplots
m = 3;
n = 2;
plotpos = @(k, l) (k-1)*n + l;
useplot = @(k, l) subplot(m, n, plotpos(k, l));

figure;

for i = 1:m
    for j = 1:n
        useplot(i, j);
        hold on;
        grid on;
    end
end

colorOrder = get(gca, 'ColorOrder');
[c, ~] = size(colorOrder);

X = zeros(N+1, M);
Y = zeros(N+1, M);

xT = zeros(1, M);
yT = zeros(1, M);

Theta = zeros(N+1, M);
A = zeros(N+1, M);

I = zeros(N+1, M);

k = 0;
for q0 = 0
    for p0 = linspace(.1, 1.9, M)
        k = k + 1;
        color = colorOrder(mod(k, c)+1, :);
        
        % St?rmer-Verlet scheme
        [p, q] = verlet(@(q) dp(0, q), p0, q0, h, N);

        % Plot
        useplot(1, 1);
        plot(q, p, '-', 'Color', color, 'UserData', k);
        
        g = N * h;
        diff_min = 0;
        start = true;
        for i = 1:length(p)
            diff = norm([p(:, i); q(:, i)] - [p0; q0]);
            if start
                if diff >= diff_min
                    diff_min = diff;
                else
                    start = false;
                end
            else
                if diff < diff_min
                    diff_min = diff;
                    g = i * h;
                else
                    break;
                end
            end
        end
        
        I(k) = i;
        
        x = H(p, q);
        y = 0:h:h*N;
        
        useplot(2, 1);
        plot(y(:, 1:i), x(:, 1:i), '-', 'Color', color, 'UserData', k);
        
        X(:, k) = x';
        Y(:, k) = y';
        
        xT(k) = x(:, i);
        yT(k) = y(:, i);
        
        theta = 2*pi / g * y;
        a = g * abs(x + 1) / (2*pi);
        
        Theta(:, k) = theta';
        A(:, k) = a';
        
        useplot(3, 1);
        plot(theta(:, 1:i), a(:, 1:i), '-', 'Color', color, 'UserData', k);
        
    end
end

useplot(3, 2);
for k = 1:M
    color = colorOrder(mod(k, c)+1, :);
    plot(Theta(1:I(k), k), A(1:I(k), k), 'Color', color, 'UserData', k);
end

X2 = zeros(N+1, M);
Y2 = zeros(N+1, M);

useplot(2, 2);
% plot(Y, X);
for k = 1:M
    color = colorOrder(mod(k, c)+1, :);
    
    theta = Theta(:, k)';
    a = A(:, k);
    g = max([I(k); 1]) * h;
    x = 2*pi*a/g - 1;
    y = theta * g/(2*pi);
    
    X2(:, k) = x';
    Y2(:, k) = y';
    
    plot(y(1:I(k)), x(1:I(k)), '-', 'Color', color, 'UserData', k);
end

useplot(1, 2);
for k = 1:M
    color = colorOrder(mod(k, c)+1, :);
    
%     p0 = sign(X2(1, k)) * sqrt(2 * abs(X2(1, k)));
    p0 = sqrt(2 * X2(1, k) + 2);
    p = zeros(1, N+1);
    q = zeros(1, N+1);
    for i = 1:N+1
        p(i) = p0 * cos(Y2(i));
        q(i) = 2 * asin(p0/2 * sin(Y2(i)));
    end
    
    plot(q, p, '-', 'Color', color, 'UserData', k);
end

for j = 1:2
    useplot(1, j);
    set(gca, 'XTick', -pi:pi/2:pi)
    set(gca, 'XTickLabel', {'-p','-p/2','0','p/2','p'}, ...
        'fontname', 'symbol')
    xlabel('q');
    ylabel('p');
    axis([-pi pi -2 2])
end

for j = 1:2
    useplot(2, j);
    plot(yT, xT, '--k');
    xlabel('y');
    ylabel('x');
%     axis([0 12 -1 1])
end

for j = 1:2
    useplot(3, j);
    set(gca, 'XTick', 0:pi/2:2*pi)
    set(gca, 'XTickLabel', {'0','p/2','p','3p/2','2p'}, ...
        'fontname', 'symbol')
    xlabel('\theta');
    ylabel('a');
    axis([0 2*pi 0 2.4])
end

% WindowButtonMotionFcn
set(gcf, 'WindowButtonDownFcn', @highlight);
% datacursormode on

range_theta = [.3 1];
range_a = [.1 2.3];

% transform_image2('images/cat24_small.png', range_theta, range_a, ...
%     @(theta, a) [ (2*pi*a-g)/g, theta*g/(2*pi) ]);

useplot(3, 2);
[im, colorMap, alpha] = imread('images/cat24_small.png');
image(range_theta, range_a, im, 'AlphaData', alpha);
colormap(colorMap);

% range_p = [.45 1.6];
% range_q = [-.5 .5];
% h = pi/3;
% N = 8;
% 
% [X, p_range, q_range, colorMap, alpha] = ...
%     transform_image('images/cat24_small.png', range_p, range_q, ...
%         @(dp, dq, p0, q0, h, N) verlet(@(q) dp(0, q), p0, q0, h, N), ...
%         dp, dq, h, N);
% 
% useplot(1, 1);
% for i = 2:2:8;
%     image(q_range{i}, p_range{i}, X{i}, 'AlphaData', alpha{i});
% end
% colormap(colorMap);
% axis([-pi pi -2 2])

% useplot(2, 1);
% axis([0 4*pi -1 1])



