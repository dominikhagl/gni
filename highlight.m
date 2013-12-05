function highlight(src, ~)
%HIGHLIGHT Mouse action function to hightlight active plots
%   All plots with the same UserData will be highlighted.
%   If UserData is 0 no plot will be highlighted.
%
%   To set UserData for a plot use
%        plot(.., .., 'UserData', k);
%   or afterwards with
%       set(h, 'UserData', k);
%   
%   To enable highlighting for the active figure use
%       set(gcf, 'WindowButtonDownFcn', @highlight);
%
%   Copyright 2013 by Dominik Hagl

    for h = findall(gcf, 'Type', 'line')
        set(h, 'LineWidth', 1);
    end

    o = get(src, 'CurrentObject');
    k = get(o, 'UserData');
    
    if k == 0
        return
    end
    
    for h = findall(gcf, 'Type', 'line', 'UserData', k)
        set(h, 'LineWidth', 3);
    end
end

