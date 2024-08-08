% Test and visualise ant_simulation

load("./antview/world5000_gray.mat")

% [x_nest, y_nest, world_bound, route_length, step_size]
test_cases = [[6.3, 8.45, 10, 7, 0.5],
                [5.4, 4.5, 10, 10, 1],
                [2.2, 3.4, 10, 5, 0.5]];


% Create figure for all subplots
figure('Position', [100, 100, 1200, 800]);

% Run tests for each case
for i = 1:size(test_cases, 1)
    x_nest = test_cases(i, 1);
    y_nest = test_cases(i, 2);
    world_bound = test_cases(i, 3);
    route_length = test_cases(i, 4);
    step_size = test_cases(i, 5);

    try
        % Call ant_simulation function
        [route, points] = ant_simulation(x_nest, y_nest, world_bound, route_length, step_size);
        
        % Perform assertions
        assert(all(points(:) >= 0 & points(:) <= world_bound), 'Points out of world bounds');
        assert(norm(route) <= route_length, 'Route length exceeds maximum');
        
        % Visualize route
        subplot(2, 3, i);
        
        % Draw the ant world
        fill(X',Y','k');
        hold on;
        
        % Plot the route
        plot(points(:,1), points(:,2), 'b-o', 'MarkerSize', 4);
        plot(x_nest, y_nest, 'r*', 'MarkerSize', 10);
        plot(points(end,1), points(end,2), 'gs', 'MarkerSize', 8); % Mark end point
        
        axis square;
        axis off;
        title(sprintf('Test Case %d', i));
    catch ME
        fprintf('Test case %d failed: %s\n', i, ME.message);
    end
end

sgtitle('Ant Simulation: Test Cases Visualisation');