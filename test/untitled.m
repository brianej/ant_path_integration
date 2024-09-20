% checks the threshold value when the view is similar or the same

load('./antview/world5000_gray.mat');
load("ant_data/nest1.mat")


% define the network parameters
numPN = 360; numKC = 20000; numEN = 1; % architecture
C_I_PN_var =5250; % input scaling parameter
g_PN_KC = 0.25; g_KC_EN = 2.0; % initialise max synaptic conductances
interval = 50; dt = 1; numTrain = num_views;
% use pre-defined network parameters
% use the weight matrix of KC-EN that stored in 'Record'
weight_matrix_KC_EN = (Record.weight_matrix_KC_EN);

% parameters for invoking ImgGrabber
eye_height = 0.01; % [m]
resolution = 400; % [degrees/pixel]
hfov = 296; % [degrees]

EN_pool = zeros(1, num_scan_img);
EN_pool_similar = zeros(1, num_scan_img);

% checks the output value when the view is similar (the view is not in the
% same place as learned but nearby)
% checks the output value when the view is the same
% Generate test views in a circular pattern around the nest center
test_views = facing_nest_views;
num_test_views = 16; % Number of test views to generate
radius = 10; % Radius of the circular pattern (in degrees)

for i = 1:size(facing_nest_views, 1)
    % Generate test views in a circular pattern around the current learned view
    for j = 1:num_test_views
        theta = (j-1) * 2 * pi / num_test_views; % Angle in radians
        test_views(end+1, 1) = facing_nest_views(i, 1) + radius * cos(theta);
        test_views(end, 2) = facing_nest_views(i, 2) + radius * sin(theta);
        test_views(end, 3) = facing_nest_views(i, 3);
    end
end

% Compute the EN_pool values for the test views
test_EN_pool = zeros(1, size(test_views, 1));

for i = 1:size(test_views, 1)
    temp_img_1 = img_grbr(test_views(i,1), test_views(i,2), eye_height, test_views(i,3), X, Y, Z, colp, fov, resolution);
    temp_img_2 = imresize(temp_img_1, [10, 36]);
    temp_img_3 = 1-double(temp_img_2)/255;
    temp_img_4 = adapthisteq(temp_img_3);
    temp_img_5 = reshape(temp_img_4, numPN, 1);
    temp_img = temp_img_5./sqrt(sum(temp_img_5.^2));
    temp_img = temp_img*C_I_PN_var;
    input = (temp_img);
    gpu_net;
    EN_count = sum(sum(spike_time_EN));
    test_EN_pool(i) = gather(EN_count);
end

% Visualize the test EN_pool values
figure;
polarscatter(test_views(:,1) * pi/180, test_views(:,2) * pi/180, 20, test_EN_pool, 'filled');
colorbar;
title('Test EN Pool Values in Circular Pattern');
