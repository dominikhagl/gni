function res = action_angle(theta, a, delta)
%ACTION_ANGLE Transport function for action angle variables
%   Transports theta and a for time delta (default: 1)

    if nargin < 3
        delta = 1;
    end

    x = 0.755905 * (...
            sqrt(10.9396*a.^2-30.9589*a+27.2638) + ...
            3.3075*a-4.68011).^(1/3) - ...
            1.32292./(sqrt(10.9396*a.^2-30.9589*a+27.2638) + ...
            3.3075*a - 4.68011 ...
        ).^(1/3);
    T = 4*ellipk(sqrt((1+x)/2));
    t = theta*T/(2*pi) + delta;
    theta = t/T * 2*pi;
    
    res = [ theta, a ];
end

