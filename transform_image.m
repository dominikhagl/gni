function [X, coords_p, coords_q, colormap, alpha] = ...
    transform_image(filename, p_range, q_range, solver, dp, dq, h, N)
%TRANSFORM_IMAGE Transform an image using a given ODE solver
%   TRANSFORM_IMAGE(filename, [p_min p_max], [q_min q_max], solver, dp, dq,
%   h, N) transforms the image represented by the file named 'filename'
%   placed into [q_min q_max] on the x-axis and  [p_min p_max] on the
%   y-axis using a given solver. The solver must be a function accepting
%   the parameters  (dp, dq, p0, q0, h, N). The image is transformed N
%   times by applying the solver with step size h to the ODE
%   dp/dt = dp, dq/dt = dq with initial values out of the ranges
%   [p_min p_max] and [q_min q_max].
%
%   The method returns a cell containing N images. coords_p and coors_q
%   spicify the limits on the y repectively on the x axis. alpha contains
%   a mask for each image to make pixels outside of the image transparent.
%
%   Each image can be plotted with
%           image(q_range{i}, p_range{i}, X{i}, 'AlphaData', alpha{i});
%
%   ATTENTION: q MUST be on the x-axis and p on the y-axis!
%
%   Copyright 2013 by Christian Buchmayr, Dominik Hagl

    [I, colormap, alpha_org] = imread(filename);
    [m, n, d] = size(I);
%     I = I(end:-1:1, :, :);
    I = flipdim(I, 1);
    alpha_org = flipdim(alpha_org, 1);
    
    if min(size(alpha_org)) == 0
        alpha_org = ones(m, n);
    end
    
    %% Initialize variables
    P = cell(1, N);
    Q = cell(1, N);
    for k = 1:N
        P{k} = zeros(m, n);
        Q{k} = zeros(m, n);
    end
    
    %% Get p0 and q0 for all pixels
    for i = 1:m
        for j = 1:n
            % Map pixel positions to p and q values
            p0 = scaleToRange(i, p_range, [1 m]);
            q0 = scaleToRange(j, q_range, [1 n]);
            
            P{1}(i, j) = p0;
            Q{1}(i, j) = q0;
        end
    end
    
    %% Solve the ODE at each pixel
    for i = 1:m
        for j = 1:n
            for k = 1:N-1
                p0 = P{k}(i, j);
                q0 = Q{k}(i, j);
                [p, q] = solver(dp, dq, p0, q0, h, 1);
                P{k+1}(i, j) = p(2);
                Q{k+1}(i, j) = q(2);
            end
        end
    end
    
    %% Create new image
    X = cell(1, N);
    coords_p = cell(1, N);
    coords_q = cell(1, N);
    alpha = cell(1, N);
    for k = 1:N
        % Get the corners of the transformed image in the coordinate system
        p_min = min(min(P{k}));
        p_max = max(max(P{k}));
        q_min = min(min(Q{k}));
        q_max = max(max(Q{k}));
        coords_p{k} = [p_min p_max]; %[p_max p_min];
        coords_q{k} = [q_min q_max];
        
        % Get the size (i.e. number of pixels) of the transformed image
        size_x = ceil(scaleToRange(q_max, [1 n], [q_min q_max]));
        size_y = ceil(scaleToRange(p_max, [1 m], [p_min p_max]));
        X{k} = 255 * ones(size_y, size_x, d, 'uint8');
        alpha{k} = false(size_y, size_x);
        
        % Get pixels from p and q values
        for i = 1:m
            for j = 1:n
                if alpha_org(i, j) == 0
                    continue
                end
                
                p = P{k}(i, j);
                q = Q{k}(i, j);
                %n+1-
                x = round(scaleToRange(q, [1 n], [q_min q_max]));
                y = round(scaleToRange(p, [1 m], [p_min p_max]));
                X{k}(y, x, :) = I(i, j, :); % get corresponding value of a pixel
                alpha{k}(y, x) = 1; % only make transformed pixels visible
            end
        end
        
        average_filter = fspecial('average', [5 5]);
        X{k} = imfilter(X{k}, average_filter, 'replicate');
        alpha{k} = ordfilt2(alpha{k}, 12, true(5));
    end
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

