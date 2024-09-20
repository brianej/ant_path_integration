% checks the threshold value when the view is similar or the same


load('world5000_gray.mat');
load("nest1.mat")

% no reward so no learning
reward = 0;

% generate another data where its a different nest view to see the en
% output

nest_2 = around_nest_views(4, 6, 0.5, 170);


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

EN_pool = zeros(1, num_views);
EN_pool_similar = zeros(1, num_views);
EN_pool_different = zeros(1, num_views);

% checks the output value when the view is the same
for i = 1:num_views % number of all images
    temp_img_1 = ImgGrabber(facing_nest_views(i,1),facing_nest_views(i,2),...
        eye_height,facing_nest_views(i,3),X,Y,Z,colp,hfov,resolution);
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
    EN_pool(i) = gather(EN_count);
end


% checks the output value when the view is similar (the view is not in the
% same place as learned but nearby)
% checks the output value when the view is the same
for i = 1:num_views % number of all images
    temp_img_1 = ImgGrabber(facing_nest_views(i,1)+0.01,facing_nest_views(i,2),...
        eye_height,facing_nest_views(i,3),X,Y,Z,colp,hfov,resolution);
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
    EN_pool_similar(i) = gather(EN_count);
end


% checks the output value when the view is different 
for i = 1:num_views % number of all images
    temp_img_1 = ImgGrabber(nest_2(i,1),nest_2(i,2),...
        eye_height,nest_2(i,3),X,Y,Z,colp,hfov,resolution);
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
    EN_pool_different(i) = gather(EN_count);
end


figure;
hold on
plot(1:numel(EN_pool), EN_pool, '-o', "DisplayName", "Same");
yline(mean(EN_pool), '-', "DisplayName", "Same");
plot(1:numel(EN_pool_similar), EN_pool_similar, '-o', "DisplayName", "Similar");
yline(mean(EN_pool_similar), '-', "DisplayName", "Similar");
plot(1:numel(EN_pool_different), EN_pool_different, '-o', "DisplayName", "Different");
yline(mean(EN_pool_different), '-', "DisplayName", "Different");
xlabel('Image Index');
ylabel('EN Count');
title('EN Pool Values');
grid on;
legend