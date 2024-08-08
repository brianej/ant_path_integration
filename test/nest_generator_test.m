% Test and visualise nest_generator function

% Parameters
world_bounds = 10; % meters
nest_radius = 0.5; % meters
num_nests = 10; % number of nests 

% Generate nests
nests = zeros(num_nests, 2);
for i = 1:num_nests
    [nests(i,1), nests(i,2)] = nest_generator(nest_radius, world_bounds);
end

% Plotting
figure;

% Draw the ant world
fill(X',Y','k');
hold on;

% Plot world bounds
rectangle('Position', [0 0 world_bounds world_bounds], 'EdgeColor', 'b', 'LineWidth', 2);

% Plot nests
for i = 1:num_nests
    x = nests(i,1);
    y = nests(i,2);
    plot(x, y, 'r.', 'MarkerSize', 20);
    
    % Plot nest radius
    theta = linspace(0, 2*pi, 100);
    circle_x = nest_radius * cos(theta) + x;
    circle_y = nest_radius * sin(theta) + y;
    plot(circle_x, circle_y, 'r-');
end

% Set axis limits and labels
xlim([0 world_bounds]);
ylim([0 world_bounds]);
xlabel('X (meters)');
ylabel('Y (meters)');
title('Nest Locations');
grid on;
axis square;

% Display legend
legend('World Bounds', 'Nest Center', 'Nest Radius');

hold off;

% Print nest coordinates
disp('Nest Coordinates (x, y):');
disp(nests);