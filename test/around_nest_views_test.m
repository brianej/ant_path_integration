% Test the function and create visualisation

load("ant_data/nest1.mat");
load('./antview/world5000_gray.mat');

nest_2 = around_nest_views(4, 6, 0.5, 0.3, 160);

% Create figure
figure;

% Draw the ant world
fill(X',Y','k');

hold on;

% Plot nest
plot(nest_x, nest_y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10);


% Plot view points and direction vectors
for i = 1:num_views
    % Plot view point
    plot(facing_nest_views(i,1), facing_nest_views(i,2), 'ro', 'MarkerFaceColor', 'r');

    % Plot direction vector
    direction_length = 0.2;
    dx = direction_length * cosd(facing_nest_views(i,3));
    dy = direction_length * sind(facing_nest_views(i,3));
    quiver(facing_nest_views(i,1), facing_nest_views(i,2), dx, dy, 1, 'k', 'LineWidth', 1.5);

    % Plot view point
    plot(nest_2(i,1), nest_2(i,2), 'ro', 'MarkerFaceColor', 'r');

    % Plot direction vector
    direction_length = 0.2;
    dx = direction_length * cosd(nest_2(i,3));
    dy = direction_length * sind(nest_2(i,3));
    quiver(nest_2(i,1), nest_2(i,2), dx, dy, 1, 'k', 'LineWidth', 1.5);
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
    fprintf(fileID, '%.2f\t%.2f\t%.2f\n', facing_nest_views(i,1), facing_nest_views(i,2), facing_nest_views(i,3));
end
fclose(fileID);

disp('Visualisation saved as "ant_views_visualisation.png"');
disp('Data saved as "ant_views_data.txt"');