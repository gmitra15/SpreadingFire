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

% Time-related variables
dt = 1; % timestep
simLength = 20; % length of simulation
numIterations = 1 + simLength/dt;

% Grid dimensions
n = 7; % width
m = 7; % length

%% Constants %%

DIRT = 0; % Dirt cell that doesn't burn
TREE = 1; % Tree cell that is not on fire
FIRE = 2; % Tree cell that is on fire

probTree = 0.75; % initial probability a cell is a tree
probFire = 0.05; % initial probability a tree is on fire

%% Initialize the grid %%
forest = zeros(n, m);

% maybe do with matrix operations instead of nested for loops?
for i = 1:n
    for j = 1:m
        if rand < probTree % note: check to see if rand is always < 1
            if rand < probFire
                forest(i, j) = FIRE;
            else
                forest(i, j) = TREE;
            end
        else
            forest(i, j) = DIRT;
        end
    end
end

disp(forest);
