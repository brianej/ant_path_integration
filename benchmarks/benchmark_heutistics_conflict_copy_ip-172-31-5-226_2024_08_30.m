% using both a vector based and image based approach(nn) 

clear all; close all; clc; reset(gpuDevice);

% load saved variables
load('./antview/world5000_gray.mat');

% no reward so no learning
reward = 0;

% control the random seed
random_seed = randi(1000);

% apply the random seed
rng(random_seed);

% step size
step_size = 0.1;

%radius
radius1 = 0.5;
radius2 = 0.3;
radius3 = 0.2;

% nest
nest = [5.1, 1];
nest_x = 5.1;
nest_y = 1;

% number of views
num_views = 170;

% nest radius
facing_nest_views_r1 = around_nest_views(nest_x, nest_y, radius1, 80);
facing_nest_views_r2 = around_nest_views(nest_x, nest_y, radius2, 50);
facing_nest_views_r3 = around_nest_views(nest_x, nest_y, radius3, 40);

% concat all views
facing_nest_views = cat(1,facing_nest_views_r1,facing_nest_views_r2);
facing_nest_views = cat(1,facing_nest_views,facing_nest_views_r3);

% define the network parameters
numPN = 360; numKC = 20000; numEN = 1; % architecture
C_I_PN_var =5250; % input scaling parameter
g_PN_KC = 0.25; g_KC_EN = 2.0; % initialise max synaptic conductances
interval = 50; dt = 1; numTrain = num_views;
connection_generator; % generate connection matrix for PN-KC and KC-EN

% parameters for invoking ImgGrabber
eye_height = 0.01; % [m]
resolution = 400; % [degrees/pixel]
hfov = 296; % [degrees]
raw_images = struct;

% convert the connection matrices into matlab-CUDA gpuArray
null_input = gpuArray.zeros(numPN,1);
connection_PN_KC = gpuArray(connection_PN_KC);
connection_KC_EN = gpuArray(connection_KC_EN);

% setup connection and weight matrix for KC-EN
weight_matrix_KC_EN = connection_KC_EN*g_KC_EN;
training_inputs = zeros(numPN, numTrain);

tStartTrain = tic;

%-------------------***     1. image generation     ***-------------------%

for i = 1:num_views
    raw_images.pos(i).origin = ImgGrabber(facing_nest_views(i,1), ...
        facing_nest_views(i,2), 0.01, facing_nest_views(i,3), X, Y, Z, colp, ...
        hfov, resolution);
end

close all;


%------------------***     2. train the network     ***-------------------%

for img = 1:numTrain % for each image in the training set
    % image pre-prosessing: reverse & adaptive histgram equalisation
    temp_img_1 = raw_images.pos(img).origin;
    temp_img_2 = imresize(temp_img_1, [10, 36]);
    temp_img_3 = 1-double(temp_img_2)/255;
    temp_img_4 = adapthisteq(temp_img_3);
    temp_img_5 = reshape(temp_img_4, numPN, 1);
    training_inputs(:,img) = temp_img_5;
end

% normalisation
training_inputs = training_inputs./(repmat(sqrt(sum(training_inputs.^2)),numPN,1));
% scale input by parameter C_I_PN_var
training_inputs = training_inputs*C_I_PN_var;
training_inputs = gpuArray(training_inputs);

clear Record; clear PN_activity KC_activity EN_activity

% initialize recordings in gpuArray
PN_activity = gpuArray.zeros(numTrain*3, 1);
KC_activity = gpuArray.zeros(numTrain*3, 1);
EN_activity = gpuArray.zeros(numTrain*3, 1);
KC_fired_total = gpuArray.zeros(numTrain*3, 1);

% Phase1 --> Pre-Training phase
for n_train_pre = 1:numTrain
    n_train_pre
    % no reward so no learning
    reward = 0;
    % call the main network file
    input = training_inputs(:,n_train_pre); % input to PN
    GPU_network_fun; % call the network
    PN_ind = find(sum(spike_time_PN, 2)>0);
    PN_count = length(PN_ind); % how many PN fired
    KC_ind = find(sum(spike_time_KC, 2)>0);
    KC_count = length(KC_ind); % how many KC fired
    KC_fired_total(KC_ind, n_train_pre) = 1;
    EN_count = sum(sum(spike_time_EN)); % how many spikes EN generated
    PN_activity(n_train_pre) = PN_count;
    KC_activity(n_train_pre) = KC_count;
    EN_activity(n_train_pre) = EN_count;
end

% Phase2 --> Training phase
for n_train = 1:numTrain
    n_train
    % reward - learning is allowed
    reward = 1;
    % call the main network file
    input = training_inputs(:,n_train); % input to PN
    GPU_network_fun; % call the network
    PN_ind = find(sum(spike_time_PN, 2)>0);
    PN_count = length(PN_ind); % how many PN fired
    KC_ind = find(sum(spike_time_KC, 2)>0);
    KC_count = length(KC_ind); % how many KC fired
    KC_fired_total(KC_ind, n_train) = 1;
    EN_count = sum(sum(spike_time_EN)); % how many spikes EN generated
    PN_activity(numTrain + n_train) = PN_count;
    KC_activity(numTrain + n_train) = KC_count;
    EN_activity(numTrain + n_train) = EN_count;
end

% Phase3 --> Post-Training phase
for n_train_post = 1:numTrain
    n_train_post
    % no reward so no learning
    reward = 0;
    % call the main network file
    input = training_inputs(:,n_train_post); % input to PN
    GPU_network_fun; % call the network
    PN_ind = find(sum(spike_time_PN, 2)>0);
    PN_count = length(PN_ind); % how many PN fired
    KC_ind = find(sum(spike_time_KC, 2)>0);
    KC_count = length(KC_ind); % how many KC fired
    KC_fired_total(KC_ind, n_train_pre) = 1;
    EN_count = sum(sum(spike_time_EN)); % how many spikes EN generated
    PN_activity(2*numTrain + n_train_post) = PN_count;
    KC_activity(2*numTrain + n_train_post) = KC_count;
    EN_activity(2*numTrain + n_train_post) = EN_count;
end

% gather information from GPU
Record.PN_activity = gather(PN_activity);
Record.KC_activity = gather(KC_activity);
Record.EN_activity = gather(EN_activity);
Record.KC_fired_total = gather(KC_fired_total);
Record.weight_matrix_KC_EN = gather(weight_matrix_KC_EN);

save("./ant_data/nest1.mat", "Record", "raw_images", "radius1", "radius2", "radius3", "num_views", "nest_x", "nest_y", "facing_nest_views");

tEndTrain = toc(tStartTrain);


%-------------------***  3. returning back to nest  ***-------------------%


% use pre-defined network parameters
% use the weight matrix of KC-EN that stored in 'Record'
weight_matrix_KC_EN = (Record.weight_matrix_KC_EN);

% distance from foraging point to nest
foraging_point = [6.3, 8.45];
nest = [5.1, 1];

count = 1;
    
% Calculate the polar coordinates
[angle_r, length_r] = cart2pol(foraging_point(1)-nest(1), foraging_point(2)-nest(2));

for i_ant = Ant1

    tStartTest = tic;

    % simulating error into the angle and 
    angle_f = (angle_r*(180/pi)-180) + randn * 5;
    
    % when ants should start using the neural model to see if the view is near the nest
    distance_near_nest = length_r - radius1 - 0.1;
    
    % Compute the direction vector components
    dx = step_size * cosd(angle_f);
    dy = step_size * sind(angle_f);
    
    % Calculate the number of steps required based on the step size
    number_steps = floor(distance_near_nest/step_size);
    
    % Generate the points using the step count
    route_home = zeros(number_steps+1, 2);  
    
    for i = 0:number_steps
        route_home(i+1, :) = foraging_point + [i * dx, i * dy];
    end
    
    
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
    
    while distance_travelled < distance_near_nest + radius1 * 2
        for i_scan = 1:num_scan_img % number of all images
            temp_heading = heading + scan_range/2 - scan_spd*(i_scan-1);
            temp_img_1 = ImgGrabber(route_home(end,1),route_home(end,2),...
                eye_height,temp_heading,X,Y,Z,colp,hfov,resolution);
            temp_img_2 = imresize(temp_img_1, [10, 36]);
            temp_img_3 = 1-double(temp_img_2)/255;
            temp_img_4 = adapthisteq(temp_img_3);
            temp_img_5 = reshape(temp_img_4, numPN, 1);
            temp_img = temp_img_5./sqrt(sum(temp_img_5.^2));
            temp_img = temp_img*C_I_PN_var;
            input = (temp_img);
    
            % call the network
            GPU_network_fun;
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
    
        true = near_nest(nest_x, nest_y, 0.1, route_home(end,1), route_home(end,2));
    
        if true
            disp "successful";
            break;
        end
    
        step_count = step_count + 1;
    end
    
    tEndTest = toc(tStartTest);

    navigation_result = struct;
    navigation_result.angle_f = angle_f;
    navigation_result.route_home = route_home;
    navigation_result.weight_matrix = Record.weight_matrix_KC_EN;
    navigation_result.Record = Record;
    navigation_result.EN_response = EN_response;
    navigation_result.distance_travelled = change_angle;
    navigation_result.tEndTest = tEndTest;
    navigation_result.tEndTrain = tendTrain;
    navigation_result.step_count = step_count;

    save(sprintf('./result/heuristics/ant%droute1.mat',count),'navigation_result');
    count = count + 1;
end

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
