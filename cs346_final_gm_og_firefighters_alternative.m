% Gautam Mitra & Owen Goldthwaite
% CS346 -- Computational Modeling and Simulation I
% May 1, 2020
%
% cs346_final_gm_og.m
% 
% Final Project: Cellular automata simulation of fire spreading
%
% To run: cs346_final_gm_og.m

%% Simulation Parameters %%%

% Seed the random number generator for testing; saves rng settings to a var
rng_set = rng(1);

% Time-related variables
dt = 1; % timestep
simLength = 200; % length of simulation
numIterations = 1 + simLength/dt;

% Grid dimensions
row_count = 100; % width
col_count = 100; % length

%% Constants %%
DRY = 1;
RAIN  = 2;

DIRT = 1; % Dirt cell that doesn't burn
GRASS = 2; % Grass cell that is not on fire
TREE = 3;  % Tree cell that is not on fire
FIRE = 4; % Cell that is on fire
WET_DIRT = 5;
WET_GRASS = 6;
WET_TREE = 7;
FIGHTER = 8;

prob_init_tree = 0.1; % initial probability a cell is a tree
prob_init_grass = 0.6; % initial probability a cell is grass
prob_init_fire = 0.025; % initial probability a tree is on fire
prob_init_fighter = 0.001; % initial probability fire fighter spawns

% number of timesteps it takes for cloud to shift one position
cloud_move_const = 1; % movement constant for rain clouds, higher is slower
rain_move_speed = 2; % Amount of cells rain can mover per timestep

% Fire Fighter constants
% Number of cells a fire fighter moves per timestep, currently must be >= 2
fighter_speed = 2; 

% Time values
initial_tree_time = 10000;
initial_grass_time = 1;
rain_wet_time = 25; % How long a cell stays wet (cant be on fire) after it rains

% Percent increase of fire to occur for each N/E/S/W tile thats on fire,
% e.g. 3 of them on fire = 3*cardinal_fire_chance_increase
cardinal_fire_chance_increase = 0.35;

% Same thing as cardinal increase but for diagonal
diag_fire_chance_increase = 0.25;

% Chance for a flaming tree to extinguish
tree_extinguish_chance = 0.1;

prob_lightning = 0.00005; % probability that a cell spontaneously ignites

% Wind constants set here, 1 is default. Currently just hard coded in so use
% realistic values, i.e if N_wind is high then S_wind should be something 
% relatively small
% Wind direction is what you would use if you were sailing / explaining the
% weather, i.e. a south wind comes from the south and therefore moves north
N_wind = 1/2;
E_wind = 1/2;
S_wind = 2;
W_wind = 1/2;

card_wind_speeds = [S_wind, W_wind, N_wind, E_wind];
diag_wind_speeds = [S_wind * W_wind, N_wind * W_wind, S_wind * E_wind, N_wind * E_wind];

% Boundary values for where to spawn fire, only spawns initially within these 
% values
fire_row_upper = 20;
fire_row_lower = 0;
fire_col_lower = 0;
fire_col_upper = 100;

% Boundary values for where to spawn rain, only spawns initially within these 
% values, should be no lower than 1 and no greater than row/col size
rain_row_lower = 1;
rain_row_upper = 25;
rain_col_lower = 35;
rain_col_upper = 65;

% List of rain boundaries that gets changed/used by rain movement
boundary_list = [rain_row_lower, rain_row_upper, rain_col_lower, rain_col_upper];

%% Set up forest grid
% Initialize forest to be all dirt
forests = ones(row_count, col_count, numIterations) * DIRT;
burn_time_grids = zeros(row_count, col_count, numIterations);
rain_grids = ones(row_count, col_count, numIterations);
wet_time_grids = zeros(row_count, col_count, numIterations);

for row = 1:row_count
    for col = 1:col_count
        % Vegetation initialization 
        if rand < prob_init_tree
            forests(row, col, 1) = TREE;
            burn_time_grids(row, col, 1) = initial_tree_time;
        elseif rand < prob_init_grass
            forests(row, col, 1) = GRASS;
            burn_time_grids(row, col, 1) = initial_grass_time;
        end

        % Fire initilization, can only spawn on vegetation
        if forests(row, col, 1) == TREE || forests(row, col, 1) == GRASS % Randomly add trees
            if rand < prob_init_fire % Set a percentage of the trees to be lit
                % Making fire spawn within fire bounds
                if( (row > fire_row_lower && row < fire_row_upper)... 
                &&  (col > fire_col_lower && row < fire_col_upper) )
                    forests(row, col, 1) = FIRE;
                end
            end
        end

        % Rain initilization
        if( (row > rain_row_lower && row < rain_row_upper)... 
        &&  (col > rain_col_lower && row < rain_col_upper) )
            rain_grids(row, col, 1) = RAIN;
        end

        % Initializing fire fighters
        if(rand < prob_init_fighter)
            forests(row, col, 1) = FIGHTER;
        end

    end
end
disp("Forest Initialized");

%% Main Simulation Loop
for frame = 2:numIterations
    isfirefighter = false;
    %% Absorbing boundary condition
    % Create a grid thats the size of the forest + 2 on each side
    extended_grid_size = size(forests( : , : , frame-1))+2;
    extended_forest = ones(extended_grid_size) * DIRT; % initialize all as dirt
    extended_rain_grid = ones(extended_grid_size) * DRY;
    
    % Set the inside portion of the grid equal to the forest values from
    % the previous timestep (a.k.a the previous frame)
    extended_forest(2:end-1, 2:end-1) = forests(:,:,frame-1);
    extended_rain_grid(2:end-1, 2:end-1) = rain_grids(:,:,frame-1);
    
    burn_time_grid = burn_time_grids(:,:,frame-1);
    wet_time_grid = wet_time_grids(:,:,frame-1);

    
    % Rain cloud movement based on the wind
    if(mod(frame, cloud_move_const) == 0 && sum([card_wind_speeds, diag_wind_speeds]) ~= 8)
        % Getting wind direction
        [val, wind_dir] = max([card_wind_speeds, diag_wind_speeds]);
        
        % Updating the rain grid
        [updated_rain_grid, b_list] = update_rain(wind_dir, extended_grid_size, boundary_list, RAIN, rain_move_speed);
        % Updating the current boundaries of the rain area
        boundary_list = b_list;
    else
        % If rain didnt update this time step
        updated_rain_grid = extended_rain_grid;
    end
    rain_grids(:, :, frame) = updated_rain_grid(2:end-1, 2:end-1);
    
    %% Loop for updating each cell in the grid
    % Loop over the indices corresponding to the original (non-extended) grid
    for row = 2:row_count + 1
        for col = 2:col_count + 1

            % Von Neumann Neighborhood
            grid_point = extended_forest(row, col);
            burn_time_point = burn_time_grid(row-1, col-1);
            wet_time_point = wet_time_grid(row-1, col-1);

            north = extended_forest(row - 1, col);
            east  = extended_forest(row, col - 1);
            south = extended_forest(row + 1, col);
            west  = extended_forest(row, col + 1);
            
            % Additional code for Moore neighborhood
            northeast = extended_forest(row - 1, col - 1);
            southeast  = extended_forest(row + 1, col - 1);
            northwest = extended_forest(row - 1, col + 1);
            southwest  = extended_forest(row + 1, col + 1);

            %% Update cell
            % If the cell was dirt, it stays dirt
            if grid_point == DIRT
                updated_grid_point = DIRT;

            % If the cell was on fire, it extinguishes in the next timestep
            elseif grid_point == FIRE
                
                % If a tile has at least 1 timestep to burn, keep it on fire
                % and decrement its remaining burn time
                if(burn_time_point > 0)
                    updated_grid_point = FIRE;
                    burn_time_point = burn_time_point - 1;
                    
                    % Chance that a TREE gets extinguished if its on fire
                    if(burn_time_point > 0 && rand < tree_extinguish_chance)
                        updated_grid_point = TREE;
                        burn_time_point = burn_time_point + 1;
                    end
                else %Tile has burned down all the way, no more burn time = dirt
                    updated_grid_point = DIRT;
                end
                
                neighbors = [north, east, south, west, northeast, southeast, northwest, southwest];

                % If neighboring a fire fighter then put turn to whatever wet
                % tile it should be. 
                if(sum(neighbors == FIGHTER) ~= 0)
                    % Wet grass if 1 burn time left
                    if(burn_time_point > 0 && burn_time_point < 2)
                        wet_time_point = rain_wet_time;
                        updated_grid_point = WET_GRASS;
                    % Wet tree if 2 or more burn time
                    elseif(burn_time_point >= 2)
                        wet_time_point = rain_wet_time;
                        updated_grid_point = WET_TREE;
                    else
                        % Otherwise dirt
                        wet_time_point = rain_wet_time;
                        updated_grid_point = WET_DIRT;
                    end
                end

            % If the cell was a tree, set it on fire or leave it alone
            elseif (grid_point == TREE || grid_point == GRASS) && wet_time_point <= 0

                % Making arrays of the direction values
                card_dirs = [north, east, south, west];
                diag_dirs = [northeast, southeast, northwest, southwest];
                
                % Logical fun to see how many neighbors are on fire then
                % multiplying number by chance of new fire
                card_fire_chance = sum((card_dirs == FIRE)...
                                        * cardinal_fire_chance_increase...
                                        .* card_wind_speeds);

                diag_fire_chance = sum((diag_dirs == FIRE)...
                                        * diag_fire_chance_increase...
                                        .* diag_wind_speeds);
                
                fire_chance = card_fire_chance + diag_fire_chance;

                if (rand < fire_chance) 
                    updated_grid_point = FIRE;
                else
                    updated_grid_point = grid_point;
                end
                    
                % Lightning strikes
                if (rand < prob_lightning)
                    updated_grid_point = FIRE;
                end    
            else
                updated_grid_point = grid_point;
            end


            rain_grid_point = rain_grids(row-1, col-1, frame);
            % Rain putting out fires and updating the wet time matrix values
            is_raining = rain_grid_point == RAIN; % if it is raining on the cell
            if(is_raining)
                wet_time_point = rain_wet_time;
                
                if(updated_grid_point == FIRE) % if on fire
                    if(grid_point == GRASS) 
                        % fire that was grass goes back to grass, essentially never
                        % igniting
                        updated_grid_point = GRASS; 
                    else
                        updated_grid_point = TREE;
                    end
                end
            else
                % Decreasing wet time if not raining on this point currently
                wet_time_point = wet_time_point - 1;
            end

            
            % Setting the cells to their wet variants if they are wet
            if(wet_time_point > 0)
                if(updated_grid_point == DIRT)
                    updated_grid_point = WET_DIRT;
                end
                
                if(updated_grid_point == GRASS)
                    updated_grid_point = WET_GRASS;
                end
                
                if(updated_grid_point == TREE)
                    updated_grid_point = WET_TREE;
                end
            else % If no longer wet return to their dry state
                if(updated_grid_point == WET_DIRT)
                    updated_grid_point = DIRT;
                end
                
                if(updated_grid_point == WET_GRASS)
                    updated_grid_point = GRASS;
                end
                
                if(updated_grid_point == WET_TREE)
                    updated_grid_point = TREE;
                end
            end

            if(grid_point == FIGHTER)
                isfirefighter = true;
                % Find the closest fire

                % Create two lists one containing row values of all tiles == to FIRE
                % and the other containing the col values
                [fire_row_locations, fire_col_locations] = ind2sub(size(forests( : , : , frame-1)), find((forests( : , : , frame-1)) == FIRE));

                % Euclidean distance formula
                % Lowest value in this list will be the "closest fire"
                row_distances = row-1 - fire_row_locations;
                col_distances = col-1 - fire_col_locations;
                distance_list = [sqrt(row_distances.^2 + col_distances.^2)];
                
                last_row = row-1;
                last_col = col-1;

                % Get the index of the closest fire to find the final position of the fire
                % we want to move too
                [val, closest_fire_idx] = min(distance_list);

                % Uses that min index to get these two row and col values from the 
                % original lists of all the fire location values
                fire_row = fire_row_locations(closest_fire_idx);
                fire_col = fire_col_locations(closest_fire_idx);


                % Variables for where the updates position of the FF will be
                new_row = row-1;
                new_col = col-1;
                
                % Moving the fire fighter
                new_row = new_row + sign(fire_row - row);
                new_col = new_col + sign(fire_col - col);


                if(new_row < 1)
                    new_row = 1;
                end

                if(new_col < 1)
                    new_col = 1;
                end
                
                updated_grid_point = DIRT;     
            end
            % Place the updated cell into the forests grid
            % Need to subtract 1 from row and 1 from col to account for the
            % forest grid being pushed down and to the right by one unit
            % from the extension
            forests(row - 1, col - 1, frame) = updated_grid_point;
            if(isfirefighter)
                forests(new_row, new_col, frame) = FIGHTER;
            end
            burn_time_grids(row - 1, col - 1, frame) = burn_time_point;
            wet_time_grids(row - 1, col - 1, frame) = wet_time_point;

        end
    end
end
disp("All forests calculated");

%% Visualize the grid

% Create the window for the animation
viz_fig = figure;
viz_axes = axes(viz_fig);

% Set the colors
dirt_color = [0.4, 0.2, 0];
grass_color = [109/255, 188/255, 0];
tree_color = [63/255, 122/255, 0];
fire_color = [237/255 41/255 57/255];
wet_dirt_color = [61/255, 47/255, 37/255];
wet_grass_color = [0/255, 188/255, 100/255];
wet_tree_color = [63/255, 122/255, 75/255];
fighter_color = [235/255, 52/255, 192/255];
figher_move_color = [0, 12/255, 255/255];

colormap(viz_axes, [dirt_color; grass_color; tree_color; fire_color; ...
wet_dirt_color; wet_grass_color; wet_tree_color; fighter_color; figher_move_color]); 

% Remove axis labels, make aspect ratio look good, and maintain that state
axis off;
axis equal;
hold on;

% rain_fig = figure;
% rain_axes = axes(rain_fig);

% rain_color = [0, 0, 1];
% dry_color  = [1, 1, 1];
% colormap(rain_axes, [dry_color; rain_color]);
% Remove axis labels, make aspect ratio look good, and maintain that state
axis off;
axis equal;
hold on;

disp("Drawing...");
for i = 1:numIterations
    w = waitforbuttonpress;
    % Turn each forest grid into an image
    image(viz_axes, forests(:, :, i));

    %image(rain_axes, rain_grids(:, :, i));

    pause(.05);
end
disp("Simulation complete!");


% Some Validation Testing Stuff

% Spawn percentages
tree_count = sum(sum(forests(:,:,1)==TREE));
fprintf('\nTree Spawn Prob: %f percent', prob_init_tree*100);
fprintf('\nTrees Spawned: %d/%d',tree_count, row_count*col_count );
fprintf('\nWhich is %f percent of cells\n', (tree_count/(row_count*col_count))*100);

grass_count = sum(sum(forests(:,:,1)==GRASS));
fprintf('\nGrass Spawn Prob: %f percent', prob_init_grass*100);
fprintf('\nGrass Spawned: %d/%d',grass_count, row_count*col_count );
fprintf('\nWhich is %f percent of cells\n', (grass_count/(row_count*col_count))*100);

fire_count = sum(sum(forests(:,:,1)==FIRE));
fprintf('NOTE: Fire can only spawn on Grass and Tree cells');
fprintf('\nFire Spawn Prob: %f percent', prob_init_fire*100);
fprintf('\nFire Spawned: %d/%d',fire_count, tree_count+grass_count );
fprintf('\nWhich is %f percent of cells\n', (fire_count/(tree_count + grass_count))*100);


function [updated_rain_grid, boundary_list] = update_rain(wind_direction, grid_size, boundary_list, RAIN, move_rate)
    updated_rain_grid = ones(grid_size);
    rain_row_lower = boundary_list(1);
    rain_row_upper = boundary_list(2);
    rain_col_lower = boundary_list(3);
    rain_col_upper = boundary_list(4);
    row_count = grid_size(1);
    col_count = grid_size(2);

    if(wind_direction == 1) % South wind, moving north
        if(rain_row_upper + move_rate <= row_count)
            rain_row_upper = rain_row_upper + move_rate;
        end
        if(rain_row_lower + move_rate <= row_count)
            rain_row_lower = rain_row_lower + move_rate;
        end
        updated_rain_grid(rain_row_lower:rain_row_upper, rain_col_lower:rain_col_upper) = RAIN;
    end

    if(wind_direction == 2) % West wind, moving east
        if(rain_col_upper + move_rate <= col_count)
            rain_col_upper = rain_col_upper + move_rate;
        end
        if(rain_col_lower + move_rate <= col_count)
            rain_col_lower = rain_col_lower + move_rate;
        end
        updated_rain_grid(rain_row_lower:rain_row_upper, rain_col_lower:rain_col_upper) = RAIN;
    end

    if(wind_direction == 3) % North wind, moving south
        if(rain_row_upper - move_rate >= 1)
            rain_row_upper = rain_row_upper - move_rate;
        end
        if(rain_row_lower - move_rate >= 1)
            rain_row_lower = rain_row_lower - move_rate;
        end
        updated_rain_grid(rain_row_lower:rain_row_upper, rain_col_lower:rain_col_upper) = RAIN;
    end
    
    if(wind_direction == 4) % East wind, moving west
        if(rain_col_upper - move_rate >= 1)
            rain_col_upper = rain_col_upper - move_rate;
        end
        if(rain_col_lower - move_rate >= 1)
            rain_col_lower = rain_col_lower - move_rate;
        end
        updated_rain_grid(rain_row_lower:rain_row_upper, rain_col_lower:rain_col_upper) = RAIN;
    end

    if(wind_direction == 5) % South West wind, moving North East
        if(rain_row_upper + move_rate <= row_count)
            rain_row_upper = rain_row_upper + move_rate;
        end
        if(rain_row_lower + move_rate <= row_count)
            rain_row_lower = rain_row_lower + move_rate;
        end
        if(rain_col_upper + move_rate <= col_count)
            rain_col_upper = rain_col_upper + move_rate;
        end
        if(rain_col_lower + move_rate <= col_count)
            rain_col_lower = rain_col_lower + move_rate;
        end
        updated_rain_grid(rain_row_lower:rain_row_upper, rain_col_lower:rain_col_upper) = RAIN;
    end

    if(wind_direction == 6) % North West wind, moving South East
        if(rain_row_upper - move_rate >= 1)
            rain_row_upper = rain_row_upper - move_rate;
        end
        if(rain_row_lower - move_rate >= 1)
            rain_row_lower = rain_row_lower - move_rate;
        end
        if(rain_col_upper + move_rate <= col_count)
            rain_col_upper = rain_col_upper + move_rate;
        end
        if(rain_col_lower + move_rate <= col_count)
            rain_col_lower = rain_col_lower + move_rate;
        end
        updated_rain_grid(rain_row_lower:rain_row_upper, rain_col_lower:rain_col_upper) = RAIN;
    end

    if(wind_direction == 7) % South East wind, moving North West
        if(rain_row_upper + move_rate <= row_count)
            rain_row_upper = rain_row_upper + move_rate;
        end
        if(rain_row_lower + move_rate <= row_count)
            rain_row_lower = rain_row_lower + move_rate;
        end
        if(rain_col_upper - move_rate >= 1)
            rain_col_upper = rain_col_upper - move_rate;
        end
        if(rain_col_lower - move_rate >= 1)
            rain_col_lower = rain_col_lower - move_rate;
        end
        updated_rain_grid(rain_row_lower:rain_row_upper, rain_col_lower:rain_col_upper) = RAIN;
    end

    if(wind_direction == 8) % North East wind, moving South West
        if(rain_row_upper - move_rate >= 1)
            rain_row_upper = rain_row_upper - move_rate;
        end
        if(rain_row_lower - move_rate >= 1)
            rain_row_lower = rain_row_lower - move_rate;
        end
        if(rain_col_upper - move_rate >= 1)
            rain_col_upper = rain_col_upper - move_rate;
        end
        if(rain_col_lower - move_rate >= 1)
            rain_col_lower = rain_col_lower - move_rate;
        end
        updated_rain_grid(rain_row_lower:rain_row_upper, rain_col_lower:rain_col_upper) = RAIN;
    end
    
    boundary_list = [rain_row_lower, rain_row_upper, rain_col_lower, rain_col_upper];
end

