close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Load the lines for visualization purposes
load('../2.0-ManualLineExtraction/lines.mat');

%% Load the metric rectification homography and variables
load('../2.2-Metric_Rectification/H_metric.mat');

%% Load the calibration matrix
load('../2.3-Intrinsic_Calibration/K.mat');

%% Load the height
load('../2.4-Height_Calculation/height.mat');

%% Compute the World position of the l lines
world_l_points = rs(:, 1:2) * l_points_metric + rs(:, 3);
world_l_points(:, 3:6) = world_l_points(:, 3:6) * l1_length/l2_length;
origin = world_l_points(:, 3);

%% Plot the S curve
load('../2.0-ManualLineExtraction/S_points.mat');
S_points = S_points_image;
S_points = [S_points; ones(1, size(S_points, 2))];
S_points = H_metric * S_points;
S_points = S_points ./ S_points(3, :);
S_points_world = rs * S_points;
scaling = 1 - height / (2 * rs(:, 3)' * r3);
S_points_world = S_points_world * scaling;
S_coordinates_ =  [rs(:, 1)' / norm(rs(:, 1)); rs(:, 2)' / norm(rs(:, 2))] * (S_points_world- origin);

%Plot the coordinates
S_color = "#8C8608";
hold on;
plot(S_coordinates_(1, :), S_coordinates_(2,:), 'Color', S_color, 'Marker', 'x', 'LineWidth', 2);
text(S_coordinates_(1, 1), S_coordinates_(2, 1) + 0.05, 'S', 'Color', S_color, 'FontSize', 12, 'FontWeight', 'bold');
plot(0, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
text(0, 0 + 0.02, 'Origin', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
axis equal;
%axis labels
xlabel('X');
ylabel('Y');
