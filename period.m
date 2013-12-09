function i = period(p)
%PERIOD Measures the period for given data
%   Detailed explanation goes here

    p0 = p(:, 1);

    diff_min = 0;
    start = true;
    
    for i = 1:length(p)
        diff = norm(p(:, i) - p0);
        if start
            if diff >= diff_min
                diff_min = diff;
            else
                start = false;
            end
        else
            if diff < diff_min
                diff_min = diff;
            else
                break;
            end
        end
    end
end

