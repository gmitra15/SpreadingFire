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
rng_set = rng(1234567);

% Time-related variables
dt = 1; % timestep
simLength = 200; % length of simulation
numIterations = 1 + simLength/dt;

frame_by_frame = false; % if turned on, you can click through animation
animation_fps = 30; % frames displayed per second in visualization

stats_mode = false; % if true, no visualization -- easier for collecting stats

% Grid dimensions
row_count = 100; % width
col_count = 100; % length

%% Constants %%
% Cell values for the rain grid
DRY = 1; % dry cell value
RAIN = 2; % raining cell value

% Cell values for the forest grid
DIRT = 1; % Dirt cell that doesn't burn
GRASS = 2; % Grass cell that is not on fire
TREE = 3;  % Tree cell that is not on fire
FIRE = 4; % Cell that is on fire
WET_DIRT = 5; % Wet dirt cell value
WET_GRASS = 6; % Wet grass cell value
WET_TREE = 7; % Wet tree cell values
FIREFIGHTER = 8; % Fire fighter cell value

%% Initialization Probabilities %%
prob_init_tree = 0.1; % initial probability a cell is a tree
prob_init_grass = 0.6; % initial probability a cell is grass
prob_init_fire = 0.0250; % initial probability a tree is on fire
prob_init_firefighter = 0.003; % initial probability fire fighter spawns

%% Fire Variables %%
% Boundary values for where to spawn fire, only spawns initially within these 
% values
fire_row_upper = 20;
fire_row_lower = 0;
fire_col_lower = 0;
fire_col_upper = 100;

% Burn duration variables
tree_burn_time = 3; % Time that a tree will burn for
grass_burn_time = 1; % Time that a grass will burn for

% Fire spreading variables
% Percent increase of fire to occur for each N/E/S/W neighbor thats on fire,
% e.g. 3 of them on fire = 3*cardinal_fire_chance_increase
cardinal_fire_chance_increase = 0.5;

% Same thing as cardinal increase but for diagonal cells (NE/NW/SE/SW)
diag_fire_chance_increase = 0.3;

% Chance for a flaming tree to extinguish on its own
tree_extinguish_prob = 0.025;

% Prob that a forest cell is struck by lightning and spotaneously ignites
prob_lightning = 0.00000;

%% Rain variables %%
% Boundary values for where to spawn rain, only spawns initially within these 
% values, should be no greater than row/col size
rain_row_lower = 35;
rain_row_upper = 65;
rain_col_lower = 1;
rain_col_upper = 15;

% List of rain boundaries that gets changed/used by rain movement
rain_bounds = [rain_row_lower, rain_row_upper, ...
                      rain_col_lower, rain_col_upper];

wet_time = 8; % How long a cell stays wet after rain or firefighter effect

% How often the rain cloud moves, higher number means slower
% e.g. rain_move_interval = 2 means the cloud will move every 2 timesteps
rain_move_interval = 1; 
% Number of cells rain can move per movement interval
rain_move_speed = 1; 

%% Wind Variables %%
% Wind constants set here, 1 is default. Currently just hard coded in so use
% realistic values, i.e if N_wind is high then S_wind should be something 
% relatively small
% Wind direction is what you would use if you were sailing / explaining the
% weather, i.e. a south wind comes from the south and therefore moves north
N_wind = 2/3;
E_wind = 2/3;
S_wind = 2/3;
W_wind = 2;

card_wind_speeds = [S_wind, W_wind, N_wind, E_wind];
diag_wind_speeds = [S_wind * W_wind, N_wind * W_wind, ...
                    S_wind * E_wind, N_wind * E_wind];

num_wind_dirs = length([card_wind_speeds, diag_wind_speeds]);

%% Set up grids %%
% Initialize forest to be all dirt
forest_grids = ones(row_count, col_count, numIterations) * DIRT;
% Intialize burn grids to be zero because nothing is burning yet
burn_time_grids = zeros(row_count, col_count, numIterations);
% Initialize rain grid to all dry
rain_grids = ones(row_count, col_count, numIterations) * DRY;
% Initialize wet grids to all zero because nothing is wet yet
wet_time_grids = zeros(row_count, col_count, numIterations);

% Loop over the grids
for row = 1:row_count
    for col = 1:col_count
    
        % Vegetation initialization
        tree_or_grass_chance = rand;
        % Setting some cells to trees
        if  tree_or_grass_chance < prob_init_tree
            forest_grids(row, col, 1) = TREE;
            burn_time_grids(row, col, 1) = tree_burn_time;
            
        % Setting some cells to grass
        % Is a sum of tree and grass prob because values from 
        % 0 to prob_init_tree become trees in the previous condition
        elseif  tree_or_grass_chance < prob_init_grass + prob_init_tree
            forest_grids(row, col, 1) = GRASS;
            burn_time_grids(row, col, 1) = grass_burn_time;
        end

        % Fire initilization, can only spawn on vegetation
        if forest_grids(row, col, 1) == TREE || ...
           forest_grids(row, col, 1) == GRASS
            if rand < prob_init_fire % Set a percentage of the veg to be lit
                % Making fire spawn within fire bounds
                if( (row > fire_row_lower && row < fire_row_upper)... 
                &&  (col > fire_col_lower && col < fire_col_upper) )
                    forest_grids(row, col, 1) = FIRE;
                end
            end
        end

        % Rain initilization
        % Making sure rain spawns within the initial boundaries set
        if( (row > rain_row_lower && row < rain_row_upper)... 
        &&  (col > rain_col_lower && col < rain_col_upper) )
            rain_grids(row, col, 1) = RAIN;
        end

        % Initializing fire fighters
        % Spawning fire fighters with their spawn prob
        if(rand < prob_init_firefighter)
            forest_grids(row, col, 1) = FIREFIGHTER;
        end

    end
end
disp("All grids initialized");

%% Main Simulation Loop
for frame = 2:numIterations
    % Boolean to check for if there is a fire fighter in this frame so no
    % random fire fighters get spawned if the workspace isnt cleared.
    % Only relevant if firefighter probability > 0 in one simulation run
    % and then is set to zero in the following run
    has_firefighter = false;
    
    %% Absorbing boundary condition
    % Create grids that are the size of the forest + 2 on each side
    extended_grid_size = size(forest_grids( : , : , frame-1))+2;
    extended_forest_grid = ones(extended_grid_size) * DIRT; % init all as dirt
    extended_rain_grid = ones(extended_grid_size) * DRY; % init all as dry
    
    % Set the inside portion of the grids equal to the corresponding values from
    % the previous timestep (a.k.a the previous frame)
    extended_forest_grid(2:end-1, 2:end-1) = forest_grids(:,:,frame-1);
    extended_rain_grid(2:end-1, 2:end-1) = rain_grids(:,:,frame-1);
    
    % Doing the same for our non-extended burn time and wet time grds
    burn_time_grid = burn_time_grids(:,:,frame-1);
    wet_time_grid = wet_time_grids(:,:,frame-1);

    
    %% Rain cloud movement based on the wind
    % If we want to move the rain clouds this frame and there is wind
    % If all cardinal wind values = 1 (i.e. no wind), this condition is skipped
    if(mod(frame, rain_move_interval) == 0 && ... 
       sum([card_wind_speeds, diag_wind_speeds]) ~= num_wind_dirs)

        % Wind direction is determined by the direction with max wind speed
        [val, wind_dir] = max([card_wind_speeds, diag_wind_speeds]);
        
        % Updating (i.e. moving) the rain grid
        [updated_rain_grid, updated_rain_bounds] = update_rain(wind_dir, ...
        extended_grid_size, rain_bounds, RAIN, rain_move_speed);
    
        % Updating the current boundaries of the rain area
        rain_bounds = updated_rain_bounds;
    else
        % If rain didnt update this time step, keep it the same
        updated_rain_grid = extended_rain_grid;
    end
    % Set new rain grid frame to the inside values (i.e. non-extended)
    rain_grids(:, :, frame) = updated_rain_grid(2:end-1, 2:end-1);
    
    %% Loop for updating each cell in the forest
    % Loop over the indices corresponding to the original (non-extended) grid
    for row = 2:row_count + 1
        for col = 2:col_count + 1

            % Store current forest grid, burn_time, and wet_time cells
            forest_cell = extended_forest_grid(row, col);

            % Subtract 1 from row, col because these grinds aren't extended
            burn_time_cell = burn_time_grid(row-1, col-1);
            wet_time_cell = wet_time_grid(row-1, col-1);
            
            % Getting Moore Neighborhood
            north = extended_forest_grid(row - 1, col);
            east  = extended_forest_grid(row, col - 1);
            south = extended_forest_grid(row + 1, col);
            west  = extended_forest_grid(row, col + 1);
            
            northeast = extended_forest_grid(row - 1, col - 1);
            southeast  = extended_forest_grid(row + 1, col - 1);
            northwest = extended_forest_grid(row - 1, col + 1);
            southwest  = extended_forest_grid(row + 1, col + 1);

            % Put all of the neighbors into a list
            neighbors = [north, east, south, west, ...
             northeast, southeast, northwest, southwest];

            %% Update cell
            % If the cell was dirt, it stays dirt
            if forest_cell == DIRT
                updated_forest_cell = DIRT;

            % If the cell was on fire, do this
            elseif forest_cell == FIRE
                
                % If a tile has at least 1 timestep to burn, keep it on fire
                % and decrement its remaining burn time
                if(burn_time_cell > 0)
                    updated_forest_cell = FIRE;
                    burn_time_cell = burn_time_cell - 1;
                    
                    % Chance that a TREE gets extinguished if its still on fire
                    if(burn_time_cell > 0 && rand < tree_extinguish_prob)
                        updated_forest_cell = TREE;
                        burn_time_cell = burn_time_cell + 1;
                    end
                else %Tile has burned down all the way, no more burn time = dirt
                    updated_forest_cell = DIRT;
                end
                
                % If neighboring a fire fighter then put out andturn to whatever 
                % the wet tile should be, based on burn time.  
                if(sum(neighbors == FIREFIGHTER) ~= 0)
                    % Set the wet time of this sell to wet_time
                    wet_time_cell = wet_time; 

                    % Determine which type of wet cell
                    % Wet grass if equal to grass burn time
                    if(burn_time_cell == grass_burn_time)
                        updated_forest_cell = WET_GRASS;
                    % Wet tree if anything longer (bc it has to be a tree)
                    elseif(burn_time_cell > grass_burn_time)
                        updated_forest_cell = WET_TREE;
                    else 
                        % Otherwise dirt
                        updated_forest_cell = WET_DIRT;
                    end
                end

            % If the cell is a tree or grass, and it isnt wet
            elseif (forest_cell == TREE || forest_cell == GRASS) ...
                    && wet_time_cell <= 0

                % Making arrays of the neighbors in groups of direction type
                card_neighbors = [north, east, south, west];
                diag_neighbors = [northeast, southeast, northwest, southwest];
                
                % Increase chance of fire by # of lit neighbors and
                % the wind directions
                card_fire_chance = sum((card_neighbors == FIRE)...
                                        * cardinal_fire_chance_increase...
                                        .* card_wind_speeds);

                diag_fire_chance = sum((diag_neighbors == FIRE)...
                                        * diag_fire_chance_increase...
                                        .* diag_wind_speeds);
                
                % sum of fire chances from cardinal and diagonal                        
                fire_chance = card_fire_chance + diag_fire_chance;

                % Light on fire if less than fire chance
                if (rand < fire_chance) 
                    updated_forest_cell = FIRE;
                else % Otherwise keep it the same as it was (tree or grass)
                    updated_forest_cell = forest_cell;
                end
                    
                % If lightning hits, then set to fire
                if (rand < prob_lightning)
                    updated_forest_cell = FIRE;
                end    
            else % Keep the same as it was
                updated_forest_cell = forest_cell;
            end

            rain_grid_cell = rain_grids(row-1, col-1, frame);
            % Rain putting out fires and updating the wet time grid values
            % if it is raining on the cell
            if(rain_grid_cell == RAIN)
                % Set the wet time value
                wet_time_cell = wet_time;

                % if cell set to fire this timestep
                if(updated_forest_cell == FIRE)
                    % If cell was grass, set to grass again
                    if(forest_cell == GRASS) 
                        updated_forest_cell = GRASS; 
                    else
                        % Otherwise, cell must be tree, so set cell to tree
                        updated_forest_cell = TREE;
                    end
                end
            else
                % Decreasing wet time if not raining on this cell currently
                wet_time_cell = wet_time_cell - 1;
            end

            
            % Setting the cells to their wet variants if they are
            % affected by rain or firefighter
            if(wet_time_cell > 0)
                if(updated_forest_cell == DIRT)
                    updated_forest_cell = WET_DIRT;
                end
                
                if(updated_forest_cell == GRASS)
                    updated_forest_cell = WET_GRASS;
                end
                
                if(updated_forest_cell == TREE)
                    updated_forest_cell = WET_TREE;
                end
            else % If no longer wet then return to their dry state
                if(updated_forest_cell == WET_DIRT)
                    updated_forest_cell = DIRT;
                end
                
                if(updated_forest_cell == WET_GRASS)
                    updated_forest_cell = GRASS;
                end
                
                if(updated_forest_cell == WET_TREE)
                    updated_forest_cell = TREE;
                end
            end

            %% If cell is a fire fighter
            if(forest_cell == FIREFIGHTER)
                % Bool to check for if there is a fire fighter this frame so no
                % random fire fighters get spawned if the workspace isnt cleared
                has_firefighter = true;

                % Find the closest fire

                % Create two lists one containing row values of all tiles
                % == to FIRE and the other containing the col values
                [fire_row_locations, fire_col_locations] = ind2sub(...
                    size(forest_grids( : , : , frame-1)),...
                    find((forest_grids( : , : , frame-1)) == FIRE));

                    
                % Calculate Euclidean distance of curent cell to all fires
                % Must adjust row, col values throughout this conditional
                % to account for grid-extension (i.e. subtract 1)
                row_distances = (row-1) - fire_row_locations;
                col_distances = (col-1) - fire_col_locations;
                % Lowest value in this list will be the "closest fire"
                distance_list = sqrt(row_distances.^2 + col_distances.^2);

                % Get the index of the closest fire to find the final position
                % of the fire we want to move to
                [val, closest_fire_idx] = min(distance_list);

                % Uses that min index to get these two row and col values from 
                % the original lists of all the fire location values
                fire_row = fire_row_locations(closest_fire_idx);
                fire_col = fire_col_locations(closest_fire_idx);

                % Variables for where the updated position of the FF will be
                % subtract 1 to account for grid-extension
                ff_dest_row = row-1;
                ff_dest_col = col-1;
                
                % Shift the FF location based on the direction from destination
                % point - current position point
                ff_dest_row = ff_dest_row + sign(fire_row - (row-1));
                ff_dest_col = ff_dest_col + sign(fire_col - (col-1));
                
                if isempty(distance_list)
                    ff_dest_row = row-1;
                    ff_dest_col = col-1;
                end

                % Constraining movement to the bounds of the forest grid
                if(ff_dest_row < 1)
                    ff_dest_row = 1;
                end

                if(ff_dest_col < 1)
                    ff_dest_col = 1;
                end

                if(ff_dest_col > col_count)
                    ff_dest_col = col_count;
                end

                if(ff_dest_row > row_count)
                    ff_dest_row = row_count;
                end
                
                % Setting cell where fire fighter was to be dirt
                updated_forest_cell = DIRT;     
            end
            % Place the updated cell into the forest_grids next grid
            % Need to subtract 1 from row and 1 from col to account for the
            % forest grid being pushed down and to the right by one unit
            % from the extension
            forest_grids(row - 1, col - 1, frame) = updated_forest_cell;
            if(has_firefighter) % If there is a fire fighter this frame
                % Update the current fire fighters position
                forest_grids(ff_dest_row, ff_dest_col, frame) = FIREFIGHTER;
            end

            % Do the same as with forest_grids for these two grids
            burn_time_grids(row - 1, col - 1, frame) = burn_time_cell;
            wet_time_grids(row - 1, col - 1, frame) = wet_time_cell;

        end
    end
end
disp("All forest_grids calculated");

if stats_mode == false
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
    fighter_color = [251/255, 255/255, 0/255];

    colormap(viz_axes, [dirt_color; grass_color; tree_color; fire_color; ...
    wet_dirt_color; wet_grass_color; wet_tree_color; fighter_color]); 

    % Remove axis labels, make aspect ratio look good, and maintain that state
    axis off;
    axis equal;
    hold on;

    disp("Drawing...");
    for i = 1:numIterations

        % Allows for frame-by-frame viewing of forest grid
        if frame_by_frame
            w = waitforbuttonpress;
        end

        title(viz_axes, 'Frame: ' + string(i))

        % Turn each forest grid into an image
        image(viz_axes, forest_grids(:, :, i));

        pause(1/animation_fps);
    end
end


disp("Simulation complete!");


% Some Validation Testing for spawn percentages based on probabilities

init_tree_count = sum(sum(forest_grids(:,:,1)==TREE));
fprintf('\nTree Spawn Prob: %f percent', prob_init_tree*100);
fprintf('\nTrees Spawned: %d/%d',init_tree_count, row_count*col_count );
fprintf('\nWhich is %f percent of cells\n', ...
        (init_tree_count/(row_count*col_count))*100);

init_grass_count = sum(sum(forest_grids(:,:,1)==GRASS));
fprintf('\nGrass Spawn Prob: %f percent', prob_init_grass*100);
fprintf('\nGrass Spawned: %d/%d',init_grass_count, row_count*col_count );
fprintf('\nWhich is %f percent of cells\n', ...
        (init_grass_count/(row_count*col_count))*100);

fire_count = sum(sum(forest_grids(:,:,1)==FIRE));
fprintf('\nNOTE: Fire can only spawn on Grass or Tree cells and within spawn');
fprintf(' boundaries. \nBased on boundaries could be significantly ');
fprintf('less than spawn prob');
fprintf('\nFire Spawn Prob: %f percent', prob_init_fire*100);
fprintf('\nFire Spawned: %d/%d',fire_count, init_tree_count+init_grass_count );
fprintf('\nWhich is %f percent of vegetation cells\n', ...
        (fire_count/(init_tree_count + init_grass_count))*100);

final_tree_count = sum(sum(forest_grids(:,:,end)==TREE));
final_grass_count = sum(sum(forest_grids(:,:,end)==GRASS));
fprintf('\n Percentage of trees burned: %.2f percent',...
    100 * (1 - final_tree_count/init_tree_count))
fprintf('\n Percentage of grass burned: %.2f percent',...
    100 * (1-final_grass_count/init_grass_count))
fprintf('\n Percentage of foliage burned: %.2f percent\n', ...
     100 * (1-((final_grass_count+final_tree_count)/...
     (init_tree_count+init_grass_count))))


% update_rain will update the boundaries of the the rain cloud based on the 
% inputted wind direction. Essentially moving the rain cloud.
% INPUT: int: wind_direction, array: grid_size (dimensions of rain grid), 
%        array: rain_bounds (current rain location), 
%        int: RAIN (rain constant), int: move_rate (rain movement rate)
%
% OUTPUT: matrix: updated_rain_grid (new locations for rain cloud),
%         array: rain_bounds (new boundaries of rain cloud)

function [updated_rain_grid, rain_bounds] = update_rain(...
          wind_direction, grid_size, rain_bounds, RAIN, move_rate)
    % Initialize new rain grid
    updated_rain_grid = ones(grid_size);
    
    % Get values of current rain boundaries
    rain_row_lower = rain_bounds(1);
    rain_row_upper = rain_bounds(2);
    rain_col_lower = rain_bounds(3);
    rain_col_upper = rain_bounds(4);
    
    % Get row and col counts
    row_count = grid_size(1);
    col_count = grid_size(2);

    % Check wind direction
    if(wind_direction == 1) % South wind, moving north
        % If movement will be valid and not out of bounds, then move in desired
        % direction for row position. Same for columns
        if(rain_row_upper + move_rate <= row_count)
            rain_row_upper = rain_row_upper + move_rate;
        end
        if(rain_row_lower + move_rate <= row_count)
            rain_row_lower = rain_row_lower + move_rate;
        end
    end

    if(wind_direction == 2) % West wind, moving east
        if(rain_col_upper + move_rate <= col_count)
            rain_col_upper = rain_col_upper + move_rate;
        end
        if(rain_col_lower + move_rate <= col_count)
            rain_col_lower = rain_col_lower + move_rate;
        end
    end

    if(wind_direction == 3) % North wind, moving south
        if(rain_row_upper - move_rate >= 1)
            rain_row_upper = rain_row_upper - move_rate;
        end
        if(rain_row_lower - move_rate >= 1)
            rain_row_lower = rain_row_lower - move_rate;
        end
    end
    
    if(wind_direction == 4) % East wind, moving west
        if(rain_col_upper - move_rate >= 1)
            rain_col_upper = rain_col_upper - move_rate;
        end
        if(rain_col_lower - move_rate >= 1)
            rain_col_lower = rain_col_lower - move_rate;
        end
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
    end

    % Update the rain grid values within the boundaries to RAIN
    updated_rain_grid(rain_row_lower:rain_row_upper, ...
                      rain_col_lower:rain_col_upper) = RAIN;    
    % Update the rain boundaries to the new values
    rain_bounds = [rain_row_lower, rain_row_upper, ...
                  rain_col_lower, rain_col_upper];
end
