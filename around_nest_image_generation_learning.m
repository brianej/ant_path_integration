% Trains the ants on snapshots around the nest with a radius of X

% load world data
load('./antview/world5000_gray.mat');

% control the random seed
random_seed = randi(1000000);
% apply the random seed
rng(random_seed);

% parameters for invoking ImgGrabber
eye_height = 0.01; % [m]
resolution = 400; % [degrees/pixel]
hfov = 296; % [degrees]
raw_images = struct;

% world bounds
world_bound = 10;

% set radius around the nest
radius = 0.5; %[m]

% num of views around the nest
num_views = 60;

% coordinates of the nest
[nest_x, nest_y] = nest_generator(radius, world_bound);

% generate views of around the nest from x(radius) away from the nest
% facing the nest
facing_nest_views = around_nest_views(nest_x, nest_y, radius, num_views);


%-------------------***     1. image generation     ***-------------------%
for i = 1:num_views
    raw_images.pos(i).origin = img_grbr(facing_nest_views(i,1), ...
        facing_nest_views(i,2), 0.01, facing_nest_views(i,3), X, Y, Z, colp, ...
        hfov, resolution);
end

save("./ant_data/nest1.mat", "raw_images", "radius", "num_views", "nest_x", "nest_y", "facing_nest_views");

close all;


%------------------***     2. train the network     ***-------------------%

% define the network parameters
numPN = 360; numKC = 20000; numEN = 1; % architecture
C_I_PN_var =5250; % input scaling parameter
g_PN_KC = 0.25; g_KC_EN = 2.0; % initialise max synaptic conductances
interval = 50; dt = 1; numTrain = num_views;
connection_generator; % generate connection matrix for PN-KC and KC-EN

% convert the connection matrices into matlab-CUDA gpuArray
null_input = zeros(numPN,1);
connection_PN_KC = (connection_PN_KC);
connection_KC_EN = (connection_KC_EN);

% setup connection and weight matrix for KC-EN
weight_matrix_KC_EN = connection_KC_EN*g_KC_EN;
training_inputs = zeros(numPN, numTrain);

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
training_inputs = (training_inputs);

clear Record; clear PN_activity KC_activity EN_activity

% initialize recordings in gpuArray
PN_activity = zeros(numTrain*3, 1);
KC_activity = zeros(numTrain*3, 1);
EN_activity = zeros(numTrain*3, 1);
KC_fired_total = zeros(numTrain*3, 1);

% Phase1 --> Pre-Training phase
for n_train_pre = 1:numTrain
    n_train_pre
    % no reward so no learning
    reward = 0;
    % call the main network file
    input = training_inputs(:,n_train_pre); % input to PN
    gpu_net; % call the network
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
    gpu_net; % call the network
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
    gpu_net; % call the network
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

save("./ant_data/nest1.mat", "Record");
