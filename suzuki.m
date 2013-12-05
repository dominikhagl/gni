function gamma = suzuki(p)
%SUZUKI Suzuki steps
%   Returns the suzuki step scaling parameters for ODE solvers of order p
%
%   Copyright 2013 by Christian Buchmayr, Dominik Hagl

    foobar = 4^( 1/(p+1) );
    gamma1 = 1 / ( 4 - foobar );
    gamma2 = foobar / ( 4 - foobar );

    gamma = [ ...
        gamma1, ...
        gamma1, ...
        -gamma2, ...
        gamma1, ...
        gamma1
    ];

end

