
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
route_length = 5;

% step size
step_size = 1;

% define the network parameters
numPN = 360; numKC = 20000; numEN = 1; % architecture
C_I_PN_var =5250; % input scaling parameter
g_PN_KC = 0.25; g_KC_EN = 2.0; % initialise max synaptic conductances
interval = 50; dt = 1; numTrain = num_views;


%-------------------***  1. generate foraging route ***-------------------%

% generate a random route for the ants
[end_vector, points] = ant_simulation(nest_x, nest_y, world_bound, route_length, step_size);

% change the end_vector into polar coordinates
[angle, length] = cart2pol(end_vector(1), end_vector(2));
angle = angle*(180/pi) - 180;

% add an error to the angle straight home
sigma = 1;
error = normrnd(0, 15);
display(error);

new_angle = angle + error;


%-------------------***  2. returning back to nest  ***-------------------%

% use pre-defined network parameters
% use the weight matrix of KC-EN that stored in 'Record'
weight_matrix_KC_EN = (Record.weight_matrix_KC_EN);

% when ants should start using the neural model to see if the view is near
% the nest
distance_near_nest = length - radius - 0.1;
x_start_search = points(end, 1) + distance_near_nest * cosd(new_angle);
y_start_search = points(end, 2) + distance_near_nest * sind(new_angle);

route_home = [points(end, 1),  points(end, 2)];
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
heading = new_angle;
distance_travelled = distance_near_nest;


% record EN output for each position
EN_pool = zeros(1, num_scan_img);

% record for all the EN responses
EN_response = struct;

while distance_travelled < length + radius * 2
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
    threshold = 12; % need test to check if this is an appropriate number
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
    step_count = step_count + 1;
end

save("./ant_data/test1.mat", "EN_response", "distance_travelled", "route_home", "end_vector", "points");


% 3. plotting the step taken on the way back on the nest %

figure;
plot(route_home(:,1), route_home(:,2), '-o');
hold on;

% Plot the nest location
plot(nest_x, nest_y, 'r*', 'MarkerSize', 10);
text(nest_x, nest_y+2, 'Nest', 'HorizontalAlignment', 'center');

% Plot a circle around the nest
circle_radius = radius;
circle_x = nest_x + circle_radius * cosd(0:360);
circle_y = nest_y + circle_radius * sind(0:360);
plot(circle_x, circle_y, '--', 'Color', [0.5 0.5 0.5]);

% Plot the initial foraging route
plot(points(:,1), points(:,2), '--', 'Color', [0.5 0.5 0.5]);

xlabel('X');
ylabel('Y');
title('Ant Route Home');
legend('Route Home', 'Nest Location', 'Nest Radius', 'Initial Foraging Route');
grid on;


% Create a figure for all subplots
figure('Position', [100, 100, 1200, 800]);

% Loop through each step starting from the first step
for i = 2:size(EN_response.activity, 2)
    % Determine the subplot position (2 rows, 3 columns)
    subplot(2, 3, i-1); % Adjusted to i-1 because subplot indices start at 1

    % Plot the response
    plot(1:numel(EN_response.activity(i).scan), EN_response.activity(i).scan, '-o', 'DisplayName', 'Same');
    hold on;
    yline(threshold, '-');
    hold off;

    % Format the plot
    title(sprintf('Step %d', i));
end

% Add a title to the entire figure
sgtitle('EN Output Going Back to Nest');