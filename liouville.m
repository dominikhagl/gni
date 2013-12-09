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

tic;
disp('Calculation action angle coordinates...');

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
        
%         i = period([p; q])
%         g = i * h;
        
        x = H(p, q);
        y = 0:h:h*N;
        
        % http://en.wikipedia.org/wiki/Pendulum_(mathematics)#Arbitrary-amplitude_period
        %g = 4*ellipk(sin(acos(-x(1))/2));
        % http://www.wolframalpha.com/input/?i=sin%28acos%28x%29%2F2%29
        g = 4*ellipk(sqrt((1+x(1))/2));
        i = round(g / h);
        
%         I(k) = i;
        
        useplot(2, 1);
        plot(y(:, 1:i), x(:, 1:i), '-', 'Color', color, 'UserData', k);
        
        X(:, k) = x';
        Y(:, k) = y';
        
        xT(k) = x(:, i);
        yT(k) = y(:, i);
        
        theta = 2*pi / g * y;
        %a = g * abs(x + 1) / (2*pi);
        % http://www.wolframalpha.com/input/?i=integral+16*y*K%28y%29+from+0+to+sqrt%28%28H%2B1%29%2F2%29
        [~, E] = ellipke(sqrt((x+1)/2));
        a = (8/9/pi) * (...
                ((3*x)+sqrt(2*(x+1))-5).*ellipk(sqrt((1+x(1))/2)) + ...
                (sqrt(2*(x+1))+8).*E ...
            );
%         a = .35*x.^3 + 1.05*x + 1.415;

        Theta(:, k) = theta';
        A(:, k) = a';
        
        useplot(3, 1);
        plot(theta(:, 1:i), a(:, 1:i), '-', 'Color', color, 'UserData', k);
        
        useplot(3, 2);
        plot(theta(:, 1:i), a(:, 1:i), '-', 'Color', color, 'UserData', k);
        
    end
end

toc
disp('Reverting calculations...');

X2 = zeros(N+1, M);
Y2 = zeros(N+1, M);

useplot(2, 2);
% plot(Y, X);
for k = 1:M
    color = colorOrder(mod(k, c)+1, :);
    
    theta = Theta(:, k)';
    a = A(:, k);
%     g = max([I(k); 1]) * h;
%     x = 2*pi*a/g - 1;
%     y = theta * g/(2*pi);
    
    x = 0.755905 * (...
            sqrt(10.9396*a.^2-30.9589*a+27.2638) + ...
            3.3075*a-4.68011).^(1/3) - ...
            1.32292./(sqrt(10.9396*a.^2-30.9589*a+27.2638) + ...
            3.3075*a - 4.68011 ...
        ).^(1/3);
    g = 4*ellipk(sqrt((1+x(1))/2));
    y = theta * g/(2*pi);

    X2(:, k) = x';
    Y2(:, k) = y';
    
    i = round(g / h);
    
    plot(y(1:i), x(1:i), '-', 'Color', color, 'UserData', k);
end

useplot(1, 2);
for k = 1:M
    color = colorOrder(mod(k, c)+1, :);
    
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

toc

%% Plot cats
disp('Plotting cats...');

how_many_cats = 4;
fprintf('0/%d', how_many_cats);

range_theta = [.3 1];
range_a = [.1 2.3];

useplot(3, 2);
[im, colormap, alpha] = imread('images/cat24.png');
imAA = flipdim(im, 1); % image in action angle coords
alphaAA = flipdim(alpha, 1);

image(range_theta, range_a, imAA, 'AlphaData', alphaAA);

e = @(y, x) [ ...
        2 * asin(sqrt(2*x + 2)/2 * sin(y)), ...
        sqrt(2*x + 2) * cos(y) ...
    ];

for l = 1:how_many_cats
    im = imAA;
    alpha = alphaAA;

    [im, range_y, range_x, ~, alpha] = ...
        transform_image2(im, colormap, alpha, range_theta, range_a, @revert_action_angle);
    useplot(2, 2);
    image(range_y, range_x, im, 'AlphaData', alpha);
    
    [im, range_q, range_p, ~, alpha] = ...
        transform_image2(im, colormap, alpha, range_y, range_x, e);
    useplot(1, 2);
    image(range_q, range_p, im, 'AlphaData', alpha);

    fprintf('\b\b\b%d/%d', l, how_many_cats);

    if l == how_many_cats
        fprintf('\b\b\b');
        break;
    end

    [imAA, range_theta, range_a, ~, alphaAA] = ...
        transform_image2(imAA, colormap, alphaAA, range_theta, range_a, @action_angle);
    useplot(3, 2);
    image(range_theta, range_a, imAA, 'AlphaData', alphaAA);
end

toc

for l = 1:n
    useplot(l, 2);
    colormap(colormap);
end
