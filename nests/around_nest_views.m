% Generates the angles around the nest needed for ants to take snapshots of
% with 2 radius 

function views = around_nest_views(nest_x, nest_y, radius, num_views)
    % the difference in angles between each views
    angle_step_size = 360 / (num_views);
    angle = 0;
    views = zeros(num_views, 3);

    for i = 1:num_views
        views(i, 1) = nest_x + radius * cosd(angle);
        views(i, 2) = nest_y + radius * sind(angle);
        % the view is facing towards the nest, hence adding a 180 degree
        % rotation
        views(i, 3) = mod(180 + angle, 360);
        angle = angle + angle_step_size;
    end
end


