% with simulated error when the ant is foraging 

% load saved variables and world
load("ant_data/nest1.mat");
load('./antview/world5000_gray.mat');

% no reward so no learning
reward = 0;

% control the random seed
random_seed = randi(1000);

% apply the random seed
rng(random_seed);

% route length
route_length = 10;

% step size
step_size = 0.25;

% define the network parameters
numPN = 360; numKC = 20000; numEN = 1; % architecture
C_I_PN_var =5250; % input scaling parameter
g_PN_KC = 0.25; g_KC_EN = 2.0; % initialise max synaptic conductances
interval = 50; dt = 1; numTrain = num_views;


%-------------------***  1. generate foraging route ***-------------------%

% generate a random route for the ants
[home_vector_r, home_vector_f, route_r, route_f] = ant_simulation_1(nest_x, nest_y, world_bound, route_length, step_size);

% change the end_vector into polar coordinates
[angle_f, length_f] = cart2pol(home_vector_f(1), home_vector_f(2));
angle_f = angle_f*(180/pi) - 180;

[angle_r, length_r] = cart2pol(home_vector_r(1), home_vector_r(2));
angle_r = angle_r*(180/pi) - 180;

% add additional error into the return home angle
angle_f = angle_f + randn * 5;


%-------------------***  2. returning back to nest  ***-------------------%

% use pre-defined network parameters
% use the weight matrix of KC-EN that stored in 'Record'
weight_matrix_KC_EN = (Record.weight_matrix_KC_EN);

% when ants should start using the neural model to see if the view is near
% the nest
distance_near_nest = length_f - radius1 - 0.1;
x_start_search = route_f(end, 1) + distance_near_nest * cosd(angle_f);
y_start_search = route_f(end, 2) + distance_near_nest * sind(angle_f);

route_home = [route_f(end, 1),  route_f(end, 2)];
route_home = [route_home; [x_start_search, y_start_search]];
step_count = 2;

% parameters for invoking ImgGrabber
eye_height = 0.01; % [m]
resolution = 400; % [degrees/pixel]
hfov = 296; % [degrees]

scan_range = 120; %[degrees]
scan_spd = 5; % [degrees]

% number of images that needed to be test (scanned)
num_scan_img = scan_range/scan_spd + 1;

step_size = 0.1; %[m]
heading = angle_f;
distance_travelled = distance_near_nest;

% record EN output for each position
EN_pool = zeros(1, num_scan_img);

% record for all the EN responses
EN_response = struct;

while distance_travelled < length_f + radius1 * 2
    for i_scan = 1:num_scan_img % number of all images
        temp_heading = heading + scan_range/2 - scan_spd*(i_scan-1);
        temp_img_1 = img_grbr(route_home(end,1),route_home(end,2),...
            eye_height,temp_heading,X,Y,Z,colp,hfov,resolution);
        temp_img_2 = imresize(temp_img_1, [10, 36]);
        temp_img_3 = 1-double(temp_img_2)/255;
        temp_img_4 = adapthisteq(temp_img_3);
        temp_img_5 = reshape(temp_img_4, numPN, 1);
        temp_img = temp_img_5./sqrt(sum(temp_img_5.^2));
        temp_img = temp_img*C_I_PN_var;
        input = (temp_img);

        % call the network
        gpu_net;
        KC_ind = find(sum(spike_time_KC, 2)>0);
        %KC_count = length(KC_ind);
        EN_count = sum(sum(spike_time_EN));
        EN_pool(i_scan) = gather(EN_count);
    end

    EN_response.activity(step_count).scan = EN_pool;
    
    % decision-making process
    threshold = 8;
    valid_directions = find(EN_pool < threshold);
    
    % if there are view/s where the EN output is bellow the threshold head
    % in the heading of the minimum EN output, else continue on from the
    % previous direction
    if ~isempty(valid_directions)
        [~, local_index] = min(EN_pool(valid_directions));
        index = valid_directions(local_index);
        heading = heading + scan_range/2 - scan_spd*(index-1);
    end
    
    % Update position
    route_home(end+1, :) = route_home(end, :) + step_size * [cosd(heading), sind(heading)];
    
    distance_travelled = distance_travelled + step_size;

    true = near_nest(nest_x, nest_y, 0.2, route_home(end,1), route_home(end,2));

    if true
        disp "successful";
        break;
    end

    step_count = step_count + 1;
end

save("./ant_data/test1.mat", "EN_response", "distance_travelled", "route_home", "home_vector_r", "home_vector_f", "route_r", "route_f");


%-------- 3. plotting the step taken on the way back on the nest ---------%

figure;
hold on;

% Plot the home route with thicker line and larger markers
plot(route_home(:,1), route_home(:,2), '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0.2 0.6 0.8]);

% Plot the nest location with a larger marker and different color
plot(nest_x, nest_y, 'r*', 'MarkerSize', 12);

% Plot the first circle around the nest with a dashed line
circle_x1 = nest_x + radius1 * cosd(0:360);
circle_y1 = nest_y + radius1 * sind(0:360);
plot(circle_x1, circle_y1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

% Plot the second circle around the nest with a different style
circle_x2 = nest_x + radius2 * cosd(0:360);
circle_y2 = nest_y + radius2 * sind(0:360);
plot(circle_x2, circle_y2, '-.', 'Color', [0.3 0.3 0.8], 'LineWidth', 1);

% Plot the initial foraging route with a dotted line and different color
plot(route_f(:,1), route_f(:,2), ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

xlabel('X', 'FontSize', 12);
ylabel('Y', 'FontSize', 12);
title('Ant Route Home', 'FontSize', 14);
legend('Route Home', 'Nest Location', 'Nest Radius 1', 'Nest Radius 2', 'Initial Foraging Route', 'Location', 'best');
grid on;
hold off;


%--------- 4. Create a figure for all subplots showing EN outputs ---------%

figure('Position', [100, 100, 1200, 800]);

% Loop through each step starting from the first step
for i = 2:size(EN_response.activity, 2)
    % Determine the subplot position (2 rows, 3 columns)
    subplot(2, 3, i-1);

    % Plot the EN response with a different color
    plot(1:numel(EN_response.activity(i).scan), EN_response.activity(i).scan, '-o', ...
         'Color', [0.1 0.5 0.3], 'LineWidth', 2, 'MarkerSize', 6);
    hold on;

    % Add the threshold line with a label
    yline(threshold, '--', 'Threshold', 'LineWidth', 1.5, 'Color', [0.8 0.2 0.2]);

    % Format the plot
    title(sprintf('Step %d', i), 'FontSize', 12);
    xlabel('Scan Index', 'FontSize', 10);
    ylabel('EN Output', 'FontSize', 10);
    grid on;
    hold off;
end

% Add a title to the entire figure
sgtitle('EN Output Going Back to Nest', 'FontSize', 16);
