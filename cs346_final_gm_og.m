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
simLength = 101; % length of simulation
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

prob_init_tree = 0.1; % initial probability a cell is a tree
prob_init_grass = 0.6; % initial probability a cell is grass
prob_init_fire = 0.025; % initial probability a tree is on fire

% number of timesteps it takes for cloud to shift one position
cloud_move_const = 1; % movement constant for rain clouds, higher is slower
rain_move_speed = 1; % Amount of cells rain can mover per timestep

% Time values
initial_tree_time = 8;
initial_grass_time = 1;

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
N_wind = 1;
E_wind = 1;
S_wind = 1;
W_wind = 1;

card_wind_speeds = [S_wind, W_wind, N_wind, E_wind];
diag_wind_speeds = [S_wind * W_wind, N_wind * W_wind, S_wind * E_wind, N_wind * E_wind];

% Boundary values for where to spawn fire, only spawns initially within these 
% values
fire_row_upper = 30;
fire_row_lower = 0;
fire_col_lower = 0;
fire_col_upper = 100;

% Boundary values for where to spawn rain, only spawns initially within these 
% values
rain_row_lower = 1;
rain_row_upper = 35;
rain_col_lower = 1;
rain_col_upper = 100;

% List of rain boundaries that gets changed/used by rain movement
boundary_list = [rain_row_lower, rain_row_upper, rain_col_lower, rain_col_upper];

%% Set up forest grid
% Initialize forest to be all dirt
forests = ones(row_count, col_count, numIterations) * DIRT;
time_grids = zeros(row_count, col_count, numIterations);
rain_grids = ones(row_count, col_count, numIterations);

for row = 1:row_count
    for col = 1:col_count

        % Vegetation initialization 
        if rand < prob_init_tree
            forests(row, col, 1) = TREE;
            time_grids(row, col, 1) = initial_tree_time;
        elseif rand < prob_init_grass
            forests(row, col, 1) = GRASS;
            time_grids(row, col, 1) = initial_grass_time;
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
    end
end
disp("Forest Initialized");

%% Main Simulation Loop
for frame = 2:numIterations
    %% Absorbing boundary condition
    % Create a grid thats the size of the forest + 2 on each side
    extended_grid_size = size(forests( : , : , frame-1))+2;
    extended_forest = ones(extended_grid_size) * DIRT; % initialize all as dirt
    extended_rain_grid = ones(extended_grid_size) * DRY;
    
    % Set the inside portion of the grid equal to the forest values from
    % the previous timestep (a.k.a the previous frame)
    extended_forest(2:end-1, 2:end-1) = forests(:,:,frame-1);
    extended_rain_grid(2:end-1, 2:end-1) = rain_grids(:,:,frame-1);
    
    time_grid = time_grids(:,:,frame-1);
    
    % Rain cloud movement based on the wind
    if(mod(frame, cloud_move_const) == 0)
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
            time_grid_point = time_grid(row-1, col-1);

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
                if(time_grid_point > 0)
                    updated_grid_point = FIRE;
                    time_grid_point = time_grid_point - 1;
                    
                    % Chance that a TREE gets extinguished if its on fire
                    if(time_grid_point > 0 && rand < tree_extinguish_chance)
                        updated_grid_point = TREE;
                        time_grid_point = time_grid_point + 1;
                    end
                else %Tile has burned down all the way, no more burn time = dirt
                    updated_grid_point = DIRT;
                end

            % If the cell was a tree, set it on fire or leave it alone
            elseif grid_point == TREE || grid_point == GRASS

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
            end
            % Place the updated cell into the forests grid
            % Need to subtract 1 from row and 1 from col to account for the
            % forest grid being pushed down and to the right by one unit
            % from the extension

            rain_grid_point = rain_grids(row-1, col-1, frame);
            % Rain putting out fires
            is_raining = rain_grid_point == RAIN; % if it is raining on the cell
            if(updated_grid_point == FIRE && is_raining) % if on fire
                if(grid_point == GRASS) 
                    % fire that was grass goes back to grass, essentially never
                    % igniting
                    updated_grid_point = GRASS; 
                else
                    % Otherwise gets set back to being to true
                    updated_grid_point = TREE;
                end
            end
            
            forests(row - 1, col - 1, frame) = updated_grid_point;
            time_grids(row - 1, col - 1, frame) = time_grid_point;
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
grass_color = [109/255, 168/255, 74/255];
tree_color = [63/255, 122/255, 77/255];
fire_color = [237/255 41/255 57/255];
colormap(viz_axes, [dirt_color; grass_color; tree_color; fire_color]); % Dirt, Tree, Fire

% Remove axis labels, make aspect ratio look good, and maintain that state
axis off;
axis equal;
hold on;

%rain_fig = figure;
rain_axes = axes(viz_fig);

rain_color = [0, 0, 1];
dry_color  = [1, 1, 1];
colormap(rain_axes, [dry_color; rain_color]);
% Remove axis labels, make aspect ratio look good, and maintain that state
axis off;
axis equal;
hold on;

disp("Drawing...");
for i = 1:numIterations
    
    % Turn each forest grid into an image
    im1 = image(viz_axes, forests(:, :, i));
    im2 = image(rain_axes, rain_grids(:, :, i), "AlphaData",.01);

    pause(0.05);
end
disp("Simulation complete!");

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