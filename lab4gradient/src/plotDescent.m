function plotDescent(f, input_filename, x_range, y_range)
    [X,Y] = meshgrid(x_range,y_range);
    Z = f(X,Y);
    figure
    hold on
    surf(X,Y,Z)
    fileID = fopen(input_filename,'r');
    formatSpec = '%f %f';
    stepsArraySize = [2 Inf];
    steps = fscanf(fileID,formatSpec,stepsArraySize);
    fclose(fileID);
    for x = steps
        plot3(x(1), x(2), f(x(1), x(2)), 'k*');
    end
end