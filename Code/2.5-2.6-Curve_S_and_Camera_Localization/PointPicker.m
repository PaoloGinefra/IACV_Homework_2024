% The purpose of this script is to manually extract some salient points from the image. These points will be useful to build the 3d model

% Load the image
im = imread('../../Assignment/Homework Image.jpg');

% Pick the points
extractedPoints = [];

figure;
imshow(im);
title('Select the points');
hold on;
pointNumber = 1;
while true
    try
        [x, y] = ginput(1);
    catch
        break;
    end
    plot(x, y, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    %Add a label to the point
    text(x, y, [num2str(pointNumber)], 'Color', 'r', 'FontSize', 14);
    point = [x; y];
    disp(point);
    extractedPoints = [extractedPoints, point];
    pointNumber = pointNumber + 1;
end

save('./points.mat', 'extractedPoints');