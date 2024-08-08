
% load saved variables in ant_data
load("ant_data/nest1.mat");

% control the random seed
random_seed = randi(1000);

% apply the random seed
rng(random_seed);

% route length
route_length = 5;

% step size
step_size = 1;

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

% when ants should start using the neural model to see if the view is near
% the nest

distance_near_nest = length - radius;
x_start_search = points(end, 1) + distance_near_nest * cos(new_angle);
y_start_search = points(end, 2) + distance_near_nest * sin(new_angle);
route_home = [[points(end,1), points(end,2)]; [x_start_search, y_start_search] ];

while true
    

end



