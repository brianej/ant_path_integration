% generates a random route for an ants simulation inside the bound of the
% world with error added on in the route

function [home_vector_r, home_vector_f, route_r, route_f] = ant_simulation_1(x_nest, y_nest, world_bound, route_length, step_size)
    % location of the nest
    nest = [x_nest, y_nest];
    
    route_r = nest;
    route_f = nest;
    
    % the total distace travelled by the ant
    distance= 0;
    
    % the relative vector relative to the ant
    home_vector_r = [0, 0];
    home_vector_f = [0, 0];
    
    % generate random initial outward angle
    heading = randi([0, 360]);
    
    % while the ant still have some distance more to  travel
    while (distance < route_length)
        % Generate a random angle change
        angle_change = randn * 60;  % Using standard normal distribution and scaling

        % Generate error
        error = randn * 0.2;
        step_size_error = step_size + error;
    
        % Add the change to the current heading
        heading = heading + angle_change;
    
        % Ensure the new heading is within [0, 360)
        heading = mod(heading, 360);
    
        % Convert to integer
        heading = round(heading);

        x_vec = step_size * cosd(heading);
        y_vec = step_size * sind(heading);
        x_vec_f = step_size_error * cosd(heading);
        y_vec_f = step_size_error * sind(heading);
    
        % checks if the ant will be going out of bounds or not, if yes it
        % will try out a new direction
        if nest(1) + home_vector_f(1) + x_vec >= world_bound || nest(2) + home_vector_f(2) + y_vec >= world_bound
            continue;
        else
            distance = distance + step_size;
            home_vector_r(1) = home_vector_r(1) + x_vec;
            home_vector_r(2) = home_vector_r(2) + y_vec;
            home_vector_f(1) = home_vector_f(1) + x_vec_f;
            home_vector_f(2) = home_vector_f(2) + y_vec_f;
            % appends the new points into the total points reached by the
            % ants
            route_r = [route_r; nest + home_vector_r];
            route_f = [route_f; nest + home_vector_f];
        end
    end
end