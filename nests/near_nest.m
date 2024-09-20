function boolean = near_nest(nest_x, nest_y, radius, x, y)
    % outputs 1 when ant is inside the area deemed near the nest
    % the point (x,y) is inside the circle
    if (nest_x - x)^2 + (nest_y - y)^2 < radius^2
        boolean = 1;
    else 
        boolean = 0;
    end
end

