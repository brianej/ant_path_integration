% grab all 15 of the path and count how long it takes to go back to nest
% and how many steps it takes. No error is simulated as any kind of error
% would not make it possible for the ant to go back home

clear all; close all; clc;

% load saved variables
load('/home/ubuntu/MATLAB-Drive/summer24_pathIntegration/antview/AntData.mat');

count = 1;

for i_ant = [Ant1, Ant2, Ant3, Ant4, Ant5, Ant6, Ant7, Ant8, Ant9, Ant10, ...
              Ant11, Ant12, Ant13, Ant14, Ant15] 

    tStart = tic;

    step_size = 0.1;

    foraging_point = [6.3, 8.45];
    nest = [5.1, 1];
        
    % Calculate the polar coordinates
    [angle_r, length_r] = cart2pol(foraging_point(1)-nest(1), foraging_point(2)-nest(2));
    
    % change to degrees
    angle_f = (angle_r*(180/pi)-180);

    display(angle_f);
    
    % Compute the direction vector components
    dx = step_size * cosd(angle_f);
    dy = step_size * sind(angle_f);
    display(dx);
    display(dy);

    % Calculate the number of steps required based on the step size
    number_steps = floor(length_r/step_size);
    
    % Generate the points using the step count
    current_position = zeros(number_steps+1, 2);  
    
    for i = 0:number_steps
        current_position(i+1, :) = foraging_point + [i * dx, i * dy];
    end
    
    % Add the final point to ensure p2 is included
    current_position = [current_position; nest];

    tEnd = toc(tStart);

    navigation_result = struct;
    navigation_result.tEnd = tEnd;
    navigation_result.step_size = step_size;
    navigation_result.current_position = current_position;
    navigation_result.number_steps = number_steps;

    save(sprintf('/home/ubuntu/MATLAB-Drive/summer24_pathIntegration/result/pure_vector/ant%droute1.mat',count),'navigation_result');

    count = count + 1;
end
