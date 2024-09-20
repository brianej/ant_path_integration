% Clear the workspace and close all figures
clear all; clc; close all;

% Load the ant world data
load('./antview/world5000_gray.mat');

% Load the ant route data
load('./antview/AntData.mat');

% Initialise the number of ants and routes
num_ants = 15;
num_routes = 2; % 

ants = [Ant1, Ant2, Ant3, Ant4, Ant5, Ant6, Ant7, Ant8, Ant9, Ant10, Ant11, Ant12, Ant13, Ant14, Ant15];

% Loop through each ant
for ant_id = 1:5

    % Id of the ant
    ants_num = ants(ant_id);

    % Create a new figure for each ant
    figure;
    
    % Draw the ant world in the background
    fill(X', Y', 'k', 'HandleVisibility','off');
    axis square;
    axis off;
    hold on;
    
    % Loop through each route for the current ant
    for route_id = 1:2
        if route_id == 1
            temp_route = ants_num.OutwardRouteData.Route1.One_cm_control_points;
        
            % Plot the route
            scatter(temp_route(:,1)/100, temp_route(:,2)/100, '.', 'DisplayName','Route 1'); %/100 to change to meters
        else
            temp_route = ants_num.OutwardRouteData.Route2.One_cm_control_points;
        
            % Plot the route
            scatter(temp_route(:,1)/100, temp_route(:,2)/100, '.', 'DisplayName','Route 2'); %/100 to change to meters
        end
    end

    % Highlight the first point (nest)
    scatter(temp_route(1,1)/100, temp_route(1,2)/100, 50, 'r', 'filled','DisplayName', 'Nest'); % Larger red dot for the nest
    
    legend
    % Add a title to the figure
    title(sprintf('Ant %d', ant_id));
end