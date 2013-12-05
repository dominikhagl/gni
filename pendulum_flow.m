figure, hold on, grid on

%q0 = 0;
for p0 = linspace(.1, 1.9, 15)
    q = @(t) 2 * asin(p0/2 * sin(t));
    p = @(t) p0 * cos(t);
    
    T = 0:.01:2*pi;
    n = length(T);
    Q = zeros(1, n);
    P = zeros(1, n);
    for i = 1:n
        t = T(i);
        Q(i) = q(t);
        P(i) = p(t);
    end
    plot(Q, P);
end