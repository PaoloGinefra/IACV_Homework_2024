close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Load the lines for visualization purposes
load('../2.0-ManualLineExtraction/lines.mat');

%% Load the metric rectification Homography and vars
load('../2.2-Metric_Rectification/H_metric.mat');

%% Load the calibration matrix
load('../2.3-Intrinsic_Calibration/K.mat');


new_r = cross(rs(:, 2), rs(:, 1));
new_r = new_r / norm(new_r) * norm(rs(:, 1));

new_rs = [rs(:, 1), new_r, rs(:, 3)];

H = K * new_rs;

H_rect = [10, 0, 0; 0 10 0; 0, 0, 1] * inv(H);

%% warp the image
tform = projective2d(H_rect');
im_warped = imwarp(im, tform);

PICK_POINTS = false;

if(PICK_POINTS)
    % plot the image
    figure;
    imshow(im_warped);
    title('Pick the three reference lines on the rectified Image');
    % Ask for input befor start
    start = input('Press Enter to start');
    
    [reference_points, reference_lines] = drawLines('Reference', 4);
    save('./vertical_reference_points.mat', 'reference_points', 'reference_lines');
else
    load('./vertical_reference_points.mat');
    figure;
    imshow(im_warped);
    hold on;
    for i = 1:size(reference_lines, 2)
        line(reference_points(1, 2*i-1:2*i), reference_points(2, 2*i-1:2*i), 'Color', 'r', 'LineWidth', 2);
        random_offset = randi(50, 1);
        text(mean(reference_points(1, 2*i-1:2*i)), mean(reference_points(2, 2*i-1:2*i)) + 20 + random_offset, ['ref ' num2str(i)], 'Color', 'r', 'FontSize', 14);
    end
end

unit = 1 / norm(reference_points(:, 1) - reference_points(:, 2));
height = unit * norm(reference_points(:, 3) - reference_points(:, 4));
inner_height = unit * norm(reference_points(:, 5) - reference_points(:, 6));
padding = (height - inner_height) / 2;
cell_width = unit * norm(reference_points(:, 7) - reference_points(:, 8));
disp('Height');
disp(height);
disp('Inner Height');
disp(inner_height);
disp('Padding');
disp(padding);
disp('Cell Width');
disp(cell_width);

save('vertical_height.mat', 'height', 'inner_height', 'padding', 'cell_width');

function [points, lines] = drawLines(name, n_lines)
points = nan(2, 2 * n_lines);
lines = nan(3, n_lines);

for i = 1:n_lines
    figure(1);
    title(['Select the end points of the ' name ' lines']);
    segment = drawline("Color", "r");
    
    points(:, 2*i-1:2*i) = segment.Position';
    a = [segment.Position(1,:)'; 1];
    b = [segment.Position(2,:)'; 1];
    l = cross(a, b);
    l = l / norm(l);
    lines(:, i) = l;
    
    midPoint = mean(segment.Position, 1);
    text(midPoint(1), midPoint(2), [name num2str(i)], 'Color', 'r', 'FontSize', 14);
end
end