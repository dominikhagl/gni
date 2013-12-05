function [X, coords_p, coords_q, colormap, alpha] = ...
    transform_image2(filename, p_range, q_range, fun)


    [I, colormap, alpha_org] = imread(filename);
    [m, n, d] = size(I);
%     I = I(end:-1:1, :, :);
    I = flipdim(I, 1);
    alpha_org = flipdim(alpha_org, 1);
    
    if min(size(alpha_org)) == 0
        alpha_org = ones(m, n);
    end
    
    %% Solve the ODE at each pixel
    for i = 1:m
        for j = 1:n
            % Map pixel positions to p and q values
            p0 = scaleToRange(i, p_range, [1 m]);
            q0 = scaleToRange(j, q_range, [1 n]);
            R = fun(p0, q0);
            P = R(1);
            Q = R(2);
        end
    end
    
    %% Create new image
    % Get the corners of the transformed image in the coordinate system
    p_min = min(min(P));
    p_max = max(max(P));
    q_min = min(min(Q));
    q_max = max(max(Q));
    coords_p = [p_min p_max]; %[p_max p_min];
    coords_q = [q_min q_max];

    % Get the size (i.e. number of pixels) of the transformed image
    size_x = ceil(scaleToRange(q_max, [1 n], [q_min q_max]));
    size_y = ceil(scaleToRange(p_max, [1 m], [p_min p_max]));
    X = 255 * ones(size_y, size_x, d, 'uint8');
    alpha = false(size_y, size_x);

    % Get pixels from p and q values
    for i = 1:m
        for j = 1:n
            if alpha_org(i, j) == 0
                continue
            end

            p = P(i, j);
            q = Q(i, j);
            %n+1-
            x = round(scaleToRange(q, [1 n], [q_min q_max]));
            y = round(scaleToRange(p, [1 m], [p_min p_max]));
            X(y, x, :) = I(i, j, :); % get corresponding value of a pixel
            alpha(y, x) = 1; % only make transformed pixels visible
        end
    end

    average_filter = fspecial('average', [5 5]);
    X = imfilter(X, average_filter, 'replicate');
    alpha = ordfilt2(alpha, 12, true(5));
end

function scaled_value = scaleToRange(value, range, origin)
%SCALETORANGE Maps a value from one range into another
%   Maps a value from range 'origin' into the range 'range'
%   
%   Copyright 2013 by Christian Buchmayr, Dominik Hagl

    origin_mid = (origin(1) + origin(2)) / 2;
    range_mid = (range(1) + range(2)) / 2;
        
	scaled_value = (value - origin_mid) / (origin(2) - origin(1)) ...
        * (range(2) - range(1)) + range_mid;
end

