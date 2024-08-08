% Test the function and create visualisation
nest_x = 5;
nest_y = 5;
radius = 2;
num_views = 60;

views = around_nest_views(nest_x, nest_y, radius, num_views);

% Create figure
figure;

% Draw the ant world
fill(X',Y','k');

hold on;

% Plot circle representing the path
theta = linspace(0, 2*pi, 100);
x_circle = nest_x + radius * cos(theta);
y_circle = nest_y + radius * sin(theta);
plot(x_circle, y_circle, 'b-');

% Plot nest
plot(nest_x, nest_y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10);

% Plot view points and direction vectors
for i = 1:num_views
    % Plot view point
    plot(views(i,1), views(i,2), 'ro', 'MarkerFaceColor', 'r');

    % Plot direction vector
    direction_length = 1;
    dx = direction_length * cosd(views(i,3));
    dy = direction_length * sind(views(i,3));
    quiver(views(i,1), views(i,2), dx, dy, 0, 'k', 'LineWidth', 1.5);
end

% Set axis limits and labels
xlim([nest_x-radius-1, nest_x+radius+1]);
ylim([nest_y-radius-1, nest_y+radius+1]);
xlabel('X');
ylabel('Y');
title('Ant Views Around Nest');
axis equal;
grid on;

% Save figure
saveas(gcf, 'ant_views_visualisation.png');

% Save data to file
fileID = fopen('ant_views_data.txt', 'w');
fprintf(fileID, 'X\tY\tHeading\n');

for i = 1:num_views
    fprintf(fileID, '%.2f\t%.2f\t%.2f\n', views(i,1), views(i,2), views(i,3));
end
fclose(fileID);

disp('Visualisation saved as "ant_views_visualisation.png"');
disp('Data saved as "ant_views_data.txt"');