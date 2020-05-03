% Gautam Mitra & Owen Goldthwaite
% CS346 -- Computational Modeling and Simulation I
% May 1, 2020
%
% cs346_final_gm_og.m
% 
% Final Project: Cellular automata simulation of fire updateing
%
% To run: cs346_final_gm_og.m

clear
close all

%% Simulation Parameters %%%

% Seed the random number generator for testing; saves rng settings to a var
rng_set = rng(1);

% Time-related variables
dt = 1; % timestep
simLength = 35; % length of simulation
numIterations = 1 + simLength/dt;

% Grid dimensions
rows = 100; % width
cols = 100; % length

%% Constants %%
DIRT = 1; % Dirt cell that doesn't burn
TREE = 2; % Tree cell that is not on fire
FIRE = 3; % Tree cell that is on fire

prob_init_tree = 0.7; % initial probability a cell is a tree
prob_init_fire = 0.025; % initial probability a tree is on fire

% chance that a tree exposed to fire from N/S/E/W ignites
prob_catch_fire_cardinal = 0.44; 

% chance that a tree  exposed to fire on a diagonal (i.e. NE/NW/SE/SW) ignites
prob_catch_fire_diag = 0.25;

prob_lightning = 0.0005; % probability that a cell spontaneously ignites

%% Set up forest grid
% Initialize forest to be all dirt
forests = ones(rows, cols, numIterations) * DIRT;
for row = 1:rows
    for col = 1:cols
        if rand < prob_init_tree % Randomly add trees
            if rand < prob_init_fire % Set a percentage of the trees to be lit
                forests(row, col, 1) = FIRE;
            else
                forests(row, col, 1) = TREE;
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
    extended_forest = ones(extended_grid_size) * DIRT; % initialize all as dirt
    
    % Set the inside portion of the grid equal to the forest values from
    % the previous timestep (a.k.a the previous frame)
    extended_forest(2:end-1, 2:end-1) = forests(:,:,frame-1);
    
    %% Loop for updating each cell in the grid
    % Loop over the indices corresponding to the original (non-extended) grid
    for row = 2:rows + 1
        for col = 2:cols + 1
            % Von Neumann Neighborhood
            grid_point  = extended_forest(row, col);
            north = extended_forest(row - 1, col);
            east  = extended_forest(row, col + 1);
            south = extended_forest(row + 1, col);
            west  = extended_forest(row, col - 1);
            
            % Additional code for Moore neighborhood
            northeast = extended_forest(row - 1, col + 1);
            southeast  = extended_forest(row + 1, col + 1);
            northwest = extended_forest(row - 1, col - 1);
            southwest  = extended_forest(row + 1, col - 1);
            
            %% Update cell
            % If the cell was dirt, it stays dirt
            if grid_point == DIRT
                updated_grid_point = DIRT;

            % If the cell was on fire, it extinguishes in the next timestep
            elseif grid_point == FIRE
                updated_grid_point = DIRT;

            % If the cell was a tree, set it on fire or leave it alone
            elseif grid_point == TREE
                
                % If a neighboring cell was on fire, it is exposed
                % Any cell exposed to fire has a probability of catching fire
                if (north==FIRE || east==FIRE || south==FIRE || west==FIRE || ... 
                        northeast==FIRE || southeast==FIRE || ...
                        northwest==FIRE || southwest==FIRE)
                    % i think it would be cooler if each adjacent (N/E/S/W) tile on
                    % fire adds to the likelihood that tree becomes on
                    % fire, and each diagonal tile adds a smaller
                    % percentage
                    if (rand < prob_catch_fire_cardinal) 
                        updated_grid_point = FIRE;
                    else
                        updated_grid_point = TREE;
                    end
                    
                % Lightning strikes
                elseif (rand < prob_lightning)
                    updated_grid_point = FIRE;
                    
                else
                    updated_grid_point = TREE;
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
tb = axtoolbar(viz_axes,{'zoomin','zoomout','rotate','restoreview'});
btn = axtoolbarbtn(tb,'state');

disp("Drawing...");
for i = 1:numIterations
    
    % End the animation early if you press the square button to the left of the
    % magnifying glass
    stop_state = get(btn, 'Value');
    if stop_state
       close(gcf);
       break;
    end
    
    % Turn each forest grid into an image
    image(viz_axes, forests(:, :, i));
    pause(0.7);
end
disp("Simulation complete!");
