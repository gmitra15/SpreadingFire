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
simLength = 35; % length of simulation
numIterations = 1 + simLength/dt;

% Grid dimensions
row_count = 25; % width
col_count = 25; % length

%% Constants %%
DIRT = 1; % Dirt cell that doesn't burn
TREE = 2; % Tree cell that is not on fire
FIRE = 3; % Tree cell that is on fire

prob_init_tree = 0.7; % initial probability a cell is a tree
prob_init_fire = 0.025; % initial probability a tree is on fire

% Percent increase of fire to occur for each N/E/S/W tile thats on fire,
% e.g. 3 of them on fire = 3*cardinal_fire_chance_increase
cardinal_fire_chance_increase = 0.35;

% Same thing as cardinal increase but for diagonal
diag_fire_chance_increase = 0.25;

prob_lightning = 0.00005; % probability that a cell spontaneously ignites

% Wind constants set here, 1 is default. Currently just hard coded in so use
% realistic values, i.e if N_wind is high then S_wind should be something 
% relatively small
% Also for this wind a North wind blows the fire towards the north, south towards
% south, East towards east, and west towards west. Not the inverse.
N_wind = 3;
E_wind = 1/2;
S_wind = 1/2;
W_wind = 1/2;

card_wind_speeds = [N_wind, W_wind, S_wind, E_wind];
diag_wind_speeds = [N_wind * W_wind, S_wind * W_wind, N_wind * E_wind, S_wind * E_wind];

% Boundary values for where to spawn fire, only spawns initially within these 
% values
fire_row_upper = 30;
fire_row_lower = 0;
fire_col_lower = 0;
fire_col_upper = 100;

%% Set up forest grid
% Initialize forest to be all dirt
%forests = ones(row_count, col_count, numIterations) * DIRT;
forests = repmat(struct('state',DIRT,'time',0), row_count, col_count);

for row = 1:row_count
    for col = 1:col_count
        if rand < prob_init_tree % Randomly add trees
            if rand < prob_init_fire % Set a percentage of the trees to be lit
                % Making fire spawn within fire bounds
                if( (row > fire_row_lower && row < fire_row_upper)... 
                &&  (col > fire_col_lower && row < fire_col_upper) )
                    tile.state = FIRE;
                    tile.time = 1;
                    forests(row, col, 1) = tile;
                end
            else
                tile.state = TREE;
                tile.time = 2;
                forests(row, col, 1) = tile;
            end
        end
    end
end
disp("Forest Initialized");

%% Main Simulation Loop
for frame = 2:numIterations
    %% Absorbing boundary condition
    % Create a grid thats the size of the forest + 2 on each side
    extended_grid_size = size(forests( : , : , frame-1))+2;
    %extended_forest = ones(extended_grid_size) * DIRT; % initialize all as dirt
    extended_forest = repmat(struct('state',DIRT,'time',0), extended_grid_size);
    
    % Set the inside portion of the grid equal to the forest values from
    % the previous timestep (a.k.a the previous frame)
    extended_forest(2:end-1, 2:end-1) = forests(:,:,frame-1);
    
    %% Loop for updating each cell in the grid
    % Loop over the indices corresponding to the original (non-extended) grid
    for row = 2:row_count + 1
        for col = 2:col_count + 1
        
            % Von Neumann Neighborhood
            grid_point  = extended_forest(row, col);
            north = extended_forest(row - 1, col).state;
            east  = extended_forest(row, col + 1).state;
            south = extended_forest(row + 1, col).state;
            west  = extended_forest(row, col - 1).state;
            
            % Additional code for Moore neighborhood
            northeast = extended_forest(row - 1, col + 1).state;
            southeast  = extended_forest(row + 1, col + 1).state;
            northwest = extended_forest(row - 1, col - 1).state;
            southwest  = extended_forest(row + 1, col - 1).state;

            %% Update cell
            % If the cell was dirt, it stays dirt
            updated_grid_point = struct('state',DIRT,'time',0);

            if grid_point.state == DIRT
                updated_grid_point.state = DIRT;

            % If the cell was on fire, it extinguishes in the next timestep
            elseif grid_point.state == FIRE
                updated_grid_point.state = DIRT;

            % If the cell was a tree, set it on fire or leave it alone
            elseif grid_point.state == TREE

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
                    updated_grid_point.state = FIRE;
                else
                    updated_grid_point.state = TREE;
                end
                    
                % Lightning strikes
                if (rand < prob_lightning)
                    updated_grid_point.state = FIRE;
                end
            end
            % Place the updated cell into the forests grid
            % Need to subtract 1 from row and 1 from col to account for the
            % forest grid being pushed down and to the right by one unit
            % from the extension
            forests(row - 1, col - 1, frame) = updated_grid_point;
        end
    end
end
disp("All forests calculated");

%% Visualize the grid

% Create the window for the animation
viz_axes = axes;

% Set the colors
dirt_color = [0.4, 0.2, 0];
tree_color = [63/255, 122/255, 77/255];
fire_color = [237/255 41/255 57/255];
colormap([dirt_color; tree_color; fire_color]); % Dirt, Tree, Fire

% Remove axis labels, make aspect ratio look good, and maintain that state
axis off;
axis equal;
hold on;

% Add a button that lets you end the program early
% tb = axtoolbar(viz_axes,{'zoomin','zoomout','rotate','restoreview'});
% btn = axtoolbarbtn(tb,'state');

disp("Drawing...");
for i = 1:numIterations
    
    % End the animation early if you press the square button to the left of the
    % magnifying glass
    % stop_state = get(btn, 'Value');
    % if stop_state
    %    close(gcf);
    %    break;
    % end
    draw_forest = ones(row_count, col_count) * DIRT; 
    for r = 1:row_count
        for c = 1:col_count
            draw_forest(r,c) = forests(r,c,i).state;
        end
    end

    % Turn each forest grid into an image
    %image(viz_axes, forests(:, :, i));
    image(viz_axes, draw_forest);

    pause(0.4);
end
disp("Simulation complete!");
