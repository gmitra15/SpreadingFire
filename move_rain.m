dt = 1; % timestep
simLength = 7; % length of simulation
numIterations = 1 + simLength/dt;



row_count = 100;
col_count = 100;

rain_matrix = ones(row_count, col_count, numIterations);
rain_matrix(40:60, 40:60, 1) = 2;

direction = "NW";

% Rain cloud moving north
if direction == "N"

    for frame = 2:numIterations
        for row = 2:row_count
            for col = 1:col_count
                rain_cell = rain_matrix(row, col, frame - 1);
                updated_rain_cell = rain_matrix(row - 1, col, frame - 1);
                rain_matrix(row, col, frame) = updated_rain_cell;
            end
        end
    end

% Rain cloud moving south
elseif direction == "S"
    for frame = 2:numIterations
        for row = 1:row_count-1
            for col = 1:col_count
                rain_cell = rain_matrix(row, col, frame - 1);
                updated_rain_cell = rain_matrix(row + 1, col, frame - 1);
                rain_matrix(row, col, frame) = updated_rain_cell;
            end
        end
    end

% Rain cloud moving east
elseif direction == "E"
    for frame = 2:numIterations
        for row = 1:row_count
            for col = 2:col_count
                rain_cell = rain_matrix(row, col, frame - 1);
                updated_rain_cell = rain_matrix(row, col-1, frame - 1);
                rain_matrix(row, col, frame) = updated_rain_cell;
            end
        end
    end
    
% Rain cloud moving west
elseif direction == "W"
    for frame = 2:numIterations
        for row = 1:row_count
            for col = 1:col_count-1
                rain_cell = rain_matrix(row, col, frame - 1);
                updated_rain_cell = rain_matrix(row, col+1, frame - 1);
                rain_matrix(row, col, frame) = updated_rain_cell;
            end
        end
    end
    
% Rain cloud moving northeast
elseif direction == "NE"
    for frame = 2:numIterations
        for row = 2:row_count
            for col = 2:col_count
                rain_cell = rain_matrix(row, col, frame - 1);
                updated_rain_cell = rain_matrix(row - 1, col - 1, frame - 1);
                rain_matrix(row, col, frame) = updated_rain_cell;
            end
        end
    end

% Rain cloud moving northwest
elseif direction == "NW"
    for frame = 2:numIterations
        for row = 2:row_count
            for col = 1:col_count-1
                rain_cell = rain_matrix(row, col, frame - 1);
                updated_rain_cell = rain_matrix(row-1, col+1, frame - 1);
                rain_matrix(row, col, frame) = updated_rain_cell;
            end
        end
    end
    
% Rain cloud moving southeast
elseif direction == "SE"
    for frame = 2:numIterations
        for row = 1:row_count-1
            for col = 2:col_count
                rain_cell = rain_matrix(row, col, frame - 1);
                updated_rain_cell = rain_matrix(row+1, col-1, frame - 1);
                rain_matrix(row, col, frame) = updated_rain_cell;
            end
        end
    end

% Rain cloud moving southwest
elseif direction == "SW"
    for frame = 2:numIterations
        for row = 1:row_count-1
            for col = 1:col_count-1
                rain_cell = rain_matrix(row, col, frame - 1);
                updated_rain_cell = rain_matrix(row+1, col+1, frame - 1);
                rain_matrix(row, col, frame) = updated_rain_cell;
            end
        end
    end
end


disp(rain_matrix(:,:,1));

rain_fig = figure;
rain_axes = axes(rain_fig);

rain_color = [0, 0, 1];
dry_color  = [1, 1, 1];
colormap(rain_axes, [dry_color; rain_color]);
axis off;
axis equal;
hold on;

disp("Drawing...");
for i = 1:numIterations
   
    % Turn each forest grid into an image
    image(rain_axes, rain_matrix(:, :, i));
    pause(0.5);
end