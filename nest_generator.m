% Generates a random nest location based on the radius and the world bound

function [x, y] = nest_generator(radius, world_bounds)
    % Convert world_bounds to cms
    world_size_cm = world_bounds * 100;
    
    % Calculate the valid range for x and y
    min_coord = radius * 100;
    max_coord = world_size_cm - min_coord;
    
    % Generate random coordinates within the valid range
    x = randi([min_coord, max_coord])/100;
    y = randi([min_coord, max_coord])/100;

end