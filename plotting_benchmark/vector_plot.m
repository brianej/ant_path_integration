clear all; close all; clc;

% Define the number of ants
num_ants = 15;

% Define the colors for plotting (using a distinct color for each ant)
colors = lines(num_ants);

% Create a new figure


% Load and plot each ant's path
for count = 1:15
    figure;
    hold on;
    % Load the result for the current ant
    load(sprintf('./result/pure_vector/ant%droute1.mat', count), 'navigation_result');
    
    % Extract data from the loaded result
    current_position = navigation_result.current_position;
    foraging_point = [6.3, 8.45];
    nest = [5.1, 1];
    
    % Plot the path of the ant
    plot(current_position(:,1), current_position(:,2), 'o-', 'Color', colors(count,:), 'DisplayName', sprintf('Ant %d Path', count));
    
    % Plot the foraging point and the nest
    plot(foraging_point(1), foraging_point(2), 'rs', 'MarkerFaceColor', 'r', 'DisplayName', 'Foraging Point');
    plot(nest(1), nest(2), 'gs', 'MarkerFaceColor', 'g', 'DisplayName', 'Nest');

    % Add labels, legend, and grid
    xlabel('X Position');
    ylabel('Y Position');
    title('Paths of Ants Returning to Nest');
    legend('show');
    grid on;
    
    % Set axis limits (optional, adjust as needed)
    xlim([0, 10]);
    ylim([0, 10]);
  
end

