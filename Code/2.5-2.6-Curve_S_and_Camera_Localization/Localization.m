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

%% load the picked points
load('./points.mat');

%% Compute the world position of the points
extractedPoints = [extractedPoints; ones(1, size(extractedPoints, 2))];
extractedPoints = H_metric * extractedPoints;
extractedPoints = extractedPoints ./ extractedPoints(3, :);
world_points = rs * extractedPoints;

world_points(:, 3:end) = world_points(:, 3:end) * l1_length / l2_length;

%% make a 3d plot of colums of rs
figure;
plot3(rs(1, 3), rs(2, 3), rs(3, 3),  '+', 'LineWidth', 2);
grid on;

% draw two arrows starting from the point and in the direction of the first 2 cols of rs
hold on;
quiver3(rs(1, 3), rs(2, 3), rs(3, 3), rs(1, 1), rs(2, 1), rs(3, 1), 0, 'LineWidth', 2);
quiver3(rs(1, 3), rs(2, 3), rs(3, 3), rs(1, 2), rs(2, 2), rs(3, 2), 0, 'LineWidth', 2);

%Plot the plane of the span of the first 2 columns of rs
%create a grid
[X, Y] = meshgrid(0:l1_length / 100:l1_length, 0:l1_length/100:0.35 * l1_length);

grid_points = [X(:), Y(:)];
grid_points_cart = rs(:, 1:2)* grid_points' + rs(:, 3);
X = reshape(grid_points_cart(1, :), size(X));
Y = reshape(grid_points_cart(2, :), size(Y));
Z = reshape(grid_points_cart(3, :), size(X));
surf(X, Y, Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

%% Plot the lines
world_l_points = rs(:, 1:2) * l_points_metric + rs(:, 3);
world_l_points(:, 3:6) = world_l_points(:, 3:6) * l1_length/l2_length;
plot3(world_l_points(1, 1:2), world_l_points(2, 1:2), world_l_points(3, 1:2), 'r', 'LineWidth', 2);
plot3(world_l_points(1, 3:4), world_l_points(2, 3:4), world_l_points(3, 3:4), 'r', 'LineWidth', 2);
plot3(world_l_points(1, 5:6), world_l_points(2, 5:6), world_l_points(3, 5:6), 'r', 'LineWidth', 2);

world_m_points = rs(:, 1:2) * m_points_metric + rs(:, 3);
world_m_points(:, 9:end) = world_m_points(:, 9:end) * l1_length/l2_length;
load('./vertical_height.mat');
m_scaling = 1 - padding / (rs(:, 3)' * r3);
world_m_points = world_m_points * m_scaling;
plot3(world_m_points(1, 1:2), world_m_points(2, 1:2), world_m_points(3, 1:2), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 3:4), world_m_points(2, 3:4), world_m_points(3, 3:4), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 5:6), world_m_points(2, 5:6), world_m_points(3, 5:6), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 7:8), world_m_points(2, 7:8), world_m_points(3, 7:8), 'g', 'LineWidth', 2);
% axis of same scale
axis equal;

% add axis labels
xlabel('X');
ylabel('Y');
zlabel('Z');

% plot a camera at the origin
plotCamera('Location', [0, 0, -0.2], 'Orientation', eye(3), 'Size', 0.1, 'Color', 'b');

plot3(world_points(1, :), world_points(2, :), world_points(3, :), 'k.', 'MarkerSize', 20);

%% Plot the S curve
load('../2.5-2.6-Curve_S_and_Camera_Localization/S_points.mat');
S_points = S_points_image;
S_points = [S_points; ones(1, size(S_points, 2))];
S_points = H_metric * S_points;
S_points = S_points ./ S_points(3, :);
S_points_world = rs * S_points;
scaling = 1 - height / (2 * rs(:, 3)' * r3);
S_points_world = S_points_world * scaling;

plot3(S_points_world(1, :), S_points_world(2, :), S_points_world(3, :), 'b', 'LineWidth', 2);

%% Move everything to the origin
origin = world_l_points(:, 3);
i_prime = rs(:, 1) / norm(rs(:, 1));
j_prime = rs(:, 2) / norm(rs(:, 2));
k_prime = r3;

R_inv = [i_prime, j_prime, k_prime, origin;
    0, 0, 0, 1];

R = inv(R_inv);

world_l_points = R * [world_l_points; ones(1, size(world_l_points, 2))];
world_m_points = R * [world_m_points; ones(1, size(world_m_points, 2))];
S_points_world = R * [S_points_world; ones(1, size(S_points_world, 2))];
camera_pos = R * [0; 0; 0; 1];
camera_direction = R * [0; 0; 1; 1];
camera_direction = camera_direction / camera_direction(4);
camera_direction = camera_direction - camera_pos;
camera_direction = camera_direction / norm(camera_direction);
camera_direction = camera_direction(1:3);
camera_direction = camera_direction / norm(camera_direction);
%Compute camera orientation from camera direction for the plot camera
i_camera_direction = [camera_direction(2); -camera_direction(1); 0];
i_camera_direction = i_camera_direction / norm(i_camera_direction);

j_camera_direction = cross(camera_direction, i_camera_direction);
cameraOrientation = inv([i_camera_direction, j_camera_direction, camera_direction]);

world_points = [world_points; ones(1, size(world_points, 2))];
world_points = R * world_points;

figure;
plot3(world_l_points(1, 1:2), world_l_points(2, 1:2), world_l_points(3, 1:2), 'r', 'LineWidth', 2);
hold on;
plot3(world_l_points(1, 3:4), world_l_points(2, 3:4), world_l_points(3, 3:4), 'r', 'LineWidth', 2);
plot3(world_l_points(1, 5:6), world_l_points(2, 5:6), world_l_points(3, 5:6), 'r', 'LineWidth', 2);

plot3(world_m_points(1, 1:2), world_m_points(2, 1:2), world_m_points(3, 1:2), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 3:4), world_m_points(2, 3:4), world_m_points(3, 3:4), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 5:6), world_m_points(2, 5:6), world_m_points(3, 5:6), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 7:8), world_m_points(2, 7:8), world_m_points(3, 7:8), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 9:10), world_m_points(2, 9:10), world_m_points(3, 9:10), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 11:12), world_m_points(2, 11:12), world_m_points(3, 11:12), 'g', 'LineWidth', 2);

plot3(S_points_world(1, :), S_points_world(2, :), S_points_world(3, :), 'b', 'LineWidth', 2);

% plot3(world_points(1, :), world_points(2, :), world_points(3, :), 'k.', 'MarkerSize', 20);
% Plot hte world points
% axis of same scale
axis equal;

% grid on
grid on;

% add axis labels
xlabel('X');
ylabel('Y');
zlabel('Z');

% plot camera noraml
camera_direction = camera_direction * 0.2;
quiver3(camera_pos(1), camera_pos(2), camera_pos(3), camera_direction(1), camera_direction(2), camera_direction(3), 0, 'LineWidth', 2);
%plot i_camera_direction
i_camera_direction = i_camera_direction * 0.2;
quiver3(camera_pos(1), camera_pos(2), camera_pos(3), i_camera_direction(1), i_camera_direction(2), i_camera_direction(3), 0, 'LineWidth', 2);
%plot j_camera_direction
j_camera_direction = j_camera_direction * 0.2;
quiver3(camera_pos(1), camera_pos(2), camera_pos(3), j_camera_direction(1), j_camera_direction(2), j_camera_direction(3), 0, 'LineWidth', 2);
% plot a camera at the origin
cam = plotCamera('Location', camera_pos(1:3), 'Orientation', cameraOrientation, 'Size', 0.05, 'Color', 'b');

%% Draw the parellelepiped
parallelepiped_color = "#BCAD5C";
parallelepiped_alpha = 0.6;
% bottom Face
bottomFace_x = [0, 1, 1, 0];
bottomFace_y = [0, 0, depth_m, depth_m];
bottomFace_z = [0, 0, 0, 0];

patch(bottomFace_x, bottomFace_y, bottomFace_z,'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

topFace_x = bottomFace_x;
topFace_y = bottomFace_y;
topFace_z = bottomFace_z + height;

patch(topFace_x, topFace_y, topFace_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

backFace_x = [0, 1, 1, 0];
backFace_y = ones(1, 4) * (depth_m);
backFace_z = [0, 0, height, height];

patch(backFace_x, backFace_y, backFace_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

leftFace_x = [0, 0, 0, 0];
leftFace_y = [0, 0, depth_m, depth_m];
leftFace_z = [0, height, height, 0];

patch(leftFace_x, leftFace_y, leftFace_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

rightFace_x = ones(1, 4);
rightFace_y = leftFace_y;
rightFace_z = leftFace_z;

patch(rightFace_x, rightFace_y, rightFace_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

inner_bottomFace_x = [padding, 1-padding, 1-padding, padding];
inner_bottomFace_y = [0, 0, depth_m - backPlateDepth, depth_m - backPlateDepth];
inner_bottomFace_z = ones(1, 4) * padding;

patch(inner_bottomFace_x, inner_bottomFace_y, inner_bottomFace_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

inner_topFace_x = inner_bottomFace_x;
inner_topFace_y = inner_bottomFace_y;
inner_topFace_z = inner_bottomFace_z + height - 2 * padding;

patch(inner_topFace_x, inner_topFace_y, inner_topFace_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

inner_backFace_x = [padding, 1-padding, 1-padding, padding];
inner_backFace_y = ones(1, 4) * (depth_m- backPlateDepth);
inner_backFace_z = [padding, padding, height-padding, height-padding];

patch(inner_backFace_x, inner_backFace_y, inner_backFace_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

inner_leftFace_x = [padding, padding, padding, padding];
inner_leftFace_y = [0, 0, depth_m- backPlateDepth, depth_m- backPlateDepth];
inner_leftFace_z = [padding, height-padding, height-padding, padding];

patch(inner_leftFace_x, inner_leftFace_y, inner_leftFace_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

inner_rightFace_x = ones(1, 4) * (1-padding);
inner_rightFace_y = inner_leftFace_y;
inner_rightFace_z = inner_leftFace_z;

patch(inner_rightFace_x, inner_rightFace_y, inner_rightFace_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

inner_first_separator_x = ones(1, 4) * (padding + cell_width);
inner_first_separator_y = [0, 0, depth_m- backPlateDepth, depth_m- backPlateDepth];
inner_first_separator_z = [padding, height-padding, height-padding, padding];

patch(inner_first_separator_x, inner_first_separator_y, inner_first_separator_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

inner_second_separator_x = ones(1, 4) * (2* padding + cell_width);
inner_second_separator_y = inner_first_separator_y;
inner_second_separator_z = inner_first_separator_z;

patch(inner_second_separator_x, inner_second_separator_y, inner_second_separator_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

inner_third_separator_x = ones(1, 4) * (2 * padding + 2*cell_width);
inner_third_separator_y = inner_first_separator_y;
inner_third_separator_z = inner_first_separator_z;

patch(inner_third_separator_x, inner_third_separator_y, inner_third_separator_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);

inner_fourth_separator_x = ones(1, 4) * (3 * padding + 2 * cell_width);
inner_fourth_separator_y = inner_first_separator_y;
inner_fourth_separator_z = inner_first_separator_z;

patch(inner_fourth_separator_x, inner_fourth_separator_y, inner_fourth_separator_z, 'r', 'FaceColor', parallelepiped_color, 'FaceAlpha', parallelepiped_alpha);