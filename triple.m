function gamma = triple(p)
%TRIPLE Triple step
%   Returns the tiple step scaling parameters for ODE solvers of order p
%
%   Copyright 2013 by Dominik Hagl, Christian Buchmayr

    foobar = 2^( 1/(p+1) );
    gamma1 = 1 / ( 2 - foobar );
    gamma2 = foobar / ( 2 - foobar );

    gamma = [ ...
        gamma1, ...
        -gamma2, ...
        gamma1
    ];
end

