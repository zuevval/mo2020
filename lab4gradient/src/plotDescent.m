function plotDescent(f, input_filename)
    fileID = fopen(input_filename,'r');
    formatSpec = '%f %f';
    stepsArraySize = [2 Inf];
    steps = fscanf(fileID,formatSpec,stepsArraySize);
    fclose(fileID);
    x_min = min(steps(1,:)) - 0.5;
    x_max = max(steps(1,:)) + 0.5;
    y_min = min(steps(2,:)) - 0.5;
    y_max = max(steps(2,:)) + 0.5;
    x_range = x_min:0.02:x_max;
    y_range = y_min:0.02:y_max;
    [X,Y] = meshgrid(x_range,y_range);
    Z = f(X,Y);
    figure
    hold on
    surf(X,Y,Z)
    for x = steps
        plot3(x(1), x(2), f(x(1), x(2)), 'k*');
    end
end