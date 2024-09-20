% generates a random route for an ants simulation inside the bound othe
% world

function [route, points] = ant_simulation(x_nest, y_nest, world_bound, route_length, step_size)
    % location of the nest
    nest = [x_nest, y_nest];
    
    points = nest;
    
    % the total distace travelled by the ant
    distance= 0;
    
    % the relative vector relative to the ant
    relative_vector = [0, 0];
    
    % generate random initial outward angle
    heading = randi([0, 360]);
    
    % while the ant still have some distance more to  travel
    while (distance < route_length)
        % Generate a random angle change
        angle_change = randn * 100;  % Using standard normal distribution and scaling
    
        % Add the change to the current heading
        heading = heading + angle_change;
    
        % Ensure the new heading is within [0, 360)
        heading = mod(heading, 360);
    
        % Convert to integer
        heading = round(heading);

        x_vec = step_size * cosd(heading);
        y_vec = step_size * sind(heading);
    
        % checks if the ant will be going out of bounds or not, if yes it
        % will try out a new direction
        if nest(1) + relative_vector(1) + x_vec >= world_bound || nest(2) + relative_vector(2) + y_vec >= world_bound
            continue;
        else
            distance = distance + step_size;
            relative_vector(1) = relative_vector(1) + x_vec;
            relative_vector(2) = relative_vector(2) + y_vec;
            % appends the new points into the total points reached by the
            % ants
            points = [points; nest + relative_vector];
        end
    end
    
    route = relative_vector;
end