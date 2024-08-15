
% load saved variables and world
load("ant_data/nest1.mat");
load('./antview/world5000_gray.mat');

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
length = sqrt(end_vector(1)^2 + end_vector(2)^2);
angle = 180 - tan(end_vector(2)/end_vector(1));

% add an error to the angle straight home
sigma = 1;
error = normrnd(0, 1);

new_angle = angle+error;


%-------------------***  2. returning back to nest  ***-------------------%

% use pre-defined network parameters
% use the weight matrix of KC-EN that stored in 'Record'
weight_matrix_KC_EN = (Record.weight_matrix_KC_EN);

% when ants should start using the neural model to see if the view is near
% the nest
distance_near_nest = length - radius;
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
scan_spd = 10; % [degrees]

% number of images that needed to be test (scanned)
num_scan_img = scan_range/scan_spd + 1;

step_size = 0.1; %[m]
heading = new_angle;
distance_travelled = distance_near_nest;

% when the ant is give or take near the nest
near_nest_x;
near_nest_y;

% record EN output for each position
EN_pool = zeros(1, num_scan_img);

% record for all the EN responses
EN_response = struct;

while distance_travelled < length + radius * 2
    for i_scan = 1:num_scan_img % number of all images
        temp_heading = heading + scan_range/2 - scan_spd*(i_scan-1);
        temp_img_1 = img_grbr(route_home(end,1),route_home(end,2),...
            eye_height,temp_heading,X,Y,Z,colp,fov,resolution);
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
        KC_count = length(KC_ind);
        EN_count = sum(sum(spike_time_EN));
        EN_pool(i_scan) = gather(EN_count);
    end

    EN_response.activity(step_count).scan = EN_pool;
    
    % decision-making process
    threshold = 0.2; % need test to check if this is an appropriate number
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



