%% Kepler problem from HLW2006, pp 8-12, 46

% Problem from p 9
dp = @(p, q) [ ...
    -q(1) / (q(1)^2 + q(2)^2)^(3/2); ...
    -q(2) / (q(1)^2 + q(2)^2)^(3/2) ...
];
dq = @(p, q) p;

% Initial values from p 12
e = .6;
 
q0 = [ 1 - e; 0 ];
p0 = [ 0; sqrt( (1+e) / (1-e) ) ];

%% Step size and number of steps
%N = 4000; % 4000 steps (p 11)
%h = .05;  % step 0.5

N = 1e2;     % 1e4 steps
h = 7.5 / N; % Interval [0, 7.5] (p 46)

%N = 5;   % Debugging values to plot step sizes
%h = .05; % (remember to set depth to 1 or 2)

%% Solve ODE
tic

% Symplectic euler
%[p, q] = seuler(dp, dq, p0, q0, h, N);

% St?rmer-Verlet scheme
% [p, q] = verlet(@(q) dp(0, q), p0, q0, h, N, true);

% St?rmer-Verlet scheme with triple jump
% [p, q] = verlet_refined(@(q) dp(0, q), p0, q0, h, N, 5, 'triple', true);

% St?rmer-Verlet scheme with suzuki jump
[p, q] = verlet_refined(@(q) dp(0, q), p0, q0, h, N, 5, 'suzuki', true);

toc

%% Plot solution
% Plot
figure;
plot(q(1,:), q(2,:));
grid on;
axis equal;

%% Display error after solving on interval [0, 7.5]
if h == 7.5/N
    % Exact solution from p 46
    q_exact = [ ...
        -.828164402690770818204757585370; ...
        .778898095658635447081654480796 ...
    ];

    p_exact = [ ...
        -.856384715343395351524486215030; ...
        -.160552150799838435254419104102 ...
    ];

    % Display error
    abs_err_q = abs(q(:, end) - q_exact);
    abs_err_p = abs(p(:, end) - p_exact);
    fprintf('Absolute error q: %0.30g\n', abs_err_q);
    fprintf('Absolute error p: %0.30g\n', abs_err_p);
end