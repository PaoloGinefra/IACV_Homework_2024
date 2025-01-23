close all;
clc;
clear;

%% Load the image
im = imread('../Assignment/Homework Image.jpg');

%% Load the lines for visualization purposes
load('./lines.mat');

%% load the line at infinity
load('./VanishingLineHorizontalPlane.mat');

%% load conic
% load('./CircleC.mat');
load('./1-FeatureExtraction/ConicExtraction/ExtractedConic.mat');

%% Compute the homography for affine rectification
l_h_inf_prime = l_h_inf_prime / l_h_inf_prime(3);
H_aff = [.1, 0, 0; 0, .1, 0; l_h_inf_prime'];
H_aff_inv = inv(H_aff);

%% Sanity check the homography
fprintf('Sanity Check [Affine Rectification]\nline at infinity after homography:');
disp(warpLine(l_h_inf_prime, H_aff_inv)');

%% Apply the homography to the image
im_warped = warpImage(im, H_aff);

%% Compute the image of the lines
m_lines_aff = warpLine(m_lines, H_aff_inv);
m_lines_aff = m_lines_aff ./ m_lines_aff(3, :);
m_points_aff = warpPoint(m_points, H_aff);

l_lines_aff = warpLine(l_lines, H_aff_inv);
l_lines_aff = l_lines_aff ./ l_lines_aff(3, :);
l_points_aff = warpPoint(l_points, H_aff);

% Sanity check: parallel lines
disp('Sanity Check [Affine Rectification]');
disp('Parallel lines after affine rectification');
m_lines_aff_angles = rad2deg(atan2(m_lines_aff(2, :), m_lines_aff(1, :)));
l_lines_aff_angles = rad2deg(atan2(l_lines_aff(2, :), l_lines_aff(1, :)));

disp('Angles of the m lines after affine rectification in degrees');
disp(m_lines_aff_angles);

disp('Angles of the l lines after affine rectification in degrees');
disp(l_lines_aff_angles);

% Compute conic Errors
conicErrors = computeConicErrors(C, im);
[center, axes, angle] = paramsFromHomogenousConic_separated(C);

%% Apply the homography
C_aff = warpConic(C, H_aff_inv);
C_aff_par = paramsFromHomogenousConic(C_aff);

%% Compute the conic errors after the affine rectification
C_aff_errors = computeConicErrors(C_aff, im_warped);

%% convert the conic coefficient to geometric parameters
[center_aff, axes_aff, angle_aff] = paramsFromHomogenousConic_separated(C_aff);

%% Metric rectification
% Rotation
U = [cos(angle_aff), -sin(angle_aff); sin(angle_aff), cos(angle_aff);];

% rescaling the axis to make them equal
a = axes_aff(1);
b = axes_aff(2);
S = diag([1, a/b]);
K = U*S*U';
H_2 = [K, zeros(2, 1); 0, 0, 1];

H_metric = H_2 * H_aff;
H_metric_inv = inv(H_metric);

%% Apply the metric rectification
im_metric = warpImage(im, H_metric);

%% Compute the image of the conic
C_metric = warpConic(C, H_metric_inv);

%% Compute the conic errors after the metric rectification
C_metric_errors = computeConicErrors(C_metric, im_metric);
[center_metric, axes_metric, angle_metric] = paramsFromHomogenousConic_separated(C_metric);

%% Compute the images of the lines
m_lines_metric = warpLine(m_lines, H_metric_inv);
m_lines_metric = m_lines_metric ./ m_lines_metric(3, :);
m_points_metric = warpPoint(m_points, H_metric);

l_lines_metric = warpLine(l_lines, H_metric_inv);
l_lines_metric = l_lines_metric ./ l_lines_metric(3, :);
l_points_metric = warpPoint(l_points, H_metric);

%% Sanity check: othogonality of rectified lines
disp('Sanity Check [Metric Rectification]: Orthogonality of rectified lines');
disp(m_lines_metric(1:2, :)' * l_lines_metric(1:2, :));

%% Show the images
figure;
subplot(2, 3, 1);
imshow(im);
title('Original Image');

subplot(2, 3, 2);
imshow(im_warped);
title('Affine Rectification');

subplot(2, 3, 3);
imshow(im_metric);
title('Metric Rectification');

subplot(2, 3, 4);
imshow(im);
hold on;
plotLines(m_points, 'g', 'm', 2);
plotLines(l_points, 'r', 'l', 2);
contour(conicErrors, [0, 0], 'b', 'LineWidth', 2);
plotConicParams(center, axes, angle, 'b');

% plot the line at infinity
xs = 1:size(im, 2);
ys = (-l_h_inf_prime(3) - l_h_inf_prime(1) * xs) / l_h_inf_prime(2);
plot(xs, ys, 'b--', 'LineWidth', 1);

title('Original Image + lines + conic');

subplot(2, 3, 5);
imshow(im_warped);
hold on;
plotLines(m_points_aff, 'g', 'm', 2);
plotLines(l_points_aff, 'r', 'l', 2);
contour(C_aff_errors, [0, 0], 'b', 'LineWidth', 2);
plotConicParams(center_aff, axes_aff, angle_aff, 'b');
title('Affine Rectification + lines + conic');

subplot(2, 3, 6);
imshow(im_metric);
hold on;
plotLines(m_points_metric, 'g', 'm', 2);
plotLines(l_points_metric, 'r', 'l', 2);
contour(C_metric_errors, [0, 0], 'b', 'LineWidth', 2);
plotConicParams(center_metric, axes_metric, angle_metric, 'b');
title('Metric Rectification + lines + conic');

%% Compute depth m
m_lines_length = vecnorm(m_points_metric(:, 1:2:end) - m_points_metric(:, 2:2:end));
l1_legnth = vecnorm(l_points_metric(:, 1) - l_points_metric(:, 2));
average_m_length = mean(m_lines_length);
depth_m = average_m_length / l1_legnth;

disp('Depth of m');
disp(depth_m);

%% Save the matric rectification homography
save('./H_metric.mat', 'H_metric');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       USEFUL FUNCTIONS                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [warpedPoint] = warpPoint(point, H)
if(size(point, 1) == 2)
    point = [point; ones(1, size(point, 2))];
end
warpedPoint = H * point;
warpedPoint = warpedPoint ./ warpedPoint(3, :);
warpedPoint = warpedPoint(1:2, :);
end

function [warpedLine] = warpLine(line, H_inv)
warpedLine = H_inv' * line;
warpedLine = warpedLine ./ vecnorm(warpedLine);
end

function [warpedConic] = warpConic(conic, H_inv)
warpedConic = H_inv' * conic * H_inv;
end

function [imageWarped] = warpImage(image, H)
h = size(image, 1);
w = size(image, 2);

corners =  [
    0, w, w, 0;
    0, 0, h, h;
    ];

corners_warped = warpPoint(corners, H);

disp('corners_warped');
disp(corners_warped);

% plot the corners
figure;
imshow(image);
hold on;
plot(corners(1, :), corners(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
for i = 1:4
    text(corners(1, i), corners(2, i), [num2str(i) ' - (' num2str(corners(1, i)) ', ' num2str(corners(2, i)) ')'], 'Color', 'r', 'FontSize', 11);
end
plot(corners_warped(1, :), corners_warped(2, :), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
for i = 1:4
    text(corners_warped(1, i), corners_warped(2, i), [num2str(i) '-(' num2str(corners_warped(1, i)) ', ' num2str(corners_warped(2, i)) ')'], 'Color', 'b', 'FontSize', 11);
end

% Plot a bounding box for the warped corners
max_x = max(corners_warped(1, :));
min_x = min(corners_warped(1, :));
max_y = max(corners_warped(2, :));
min_y = min(corners_warped(2, :));

line([min_x, min_x], [min_y, max_y], 'Color', 'g', 'LineWidth', 2);
line([max_x, max_x], [min_y, max_y], 'Color', 'g', 'LineWidth', 2);
line([min_x, max_x], [min_y, min_y], 'Color', 'g', 'LineWidth', 2);
line([min_x, max_x], [max_y, max_y], 'Color', 'g', 'LineWidth', 2);
hold off;


h_warped = ceil(max(corners_warped(2, :)) - min(corners_warped(2, :)));
w_warped = ceil(max(corners_warped(1, :)) - min(corners_warped(1, :)));
disp('w_warped');
disp(w_warped);
disp('h_warped');
disp(h_warped);

imageWarped = zeros(h_warped, w_warped, 3, 'uint8');

for i = 1:size(imageWarped, 1)
    for j = 1:size(imageWarped, 2)
        point = [j, i, 1];
        point_original = H\(point');
        point_original = point_original ./ point_original(3);
        x = point_original(1);
        y = point_original(2);
        if x >= 1 && x <= size(image, 2) && y >= 1 && y <= size(image, 1)
            y1 = floor(y);
            y2 = ceil(y);
            
            x1 = floor(x);
            x2 = ceil(x);
            
            v1 = image(y1, x1, :) .* (x2 - x) + image(y1, x2, :) .* (x - x1);
            v2 = image(y2, x1, :) .* (x2 - x) + image(y2, x2, :) .* (x - x1);
            imageWarped(i, j, :) = v1 .* (y2 - y) + v2 .* (y - y1);
        end
    end
end
end

function [conicErrors] = computeConicErrors(C, im)
conicErrors = zeros(size(im, 1), size(im, 2));

for i = 1:size(im, 1)
    for j = 1:size(im, 2)
        point = [j, i, 1];
        conicErrors(i, j) = point * C * point';
    end
end
end

function [conicParams] = paramsFromHomogenousConic(C)
conicParams = AtoG([C(1, 1), 2 * C(1, 2), C(2, 2), 2 * C(1, 3), 2 * C(2, 3), C(3, 3)]);
end

function [center, axes, angle] = paramsFromHomogenousConic_separated(C)
conicParams = paramsFromHomogenousConic(C);
center = conicParams(1:2);
axes = conicParams(3:4);
angle = conicParams(5);
end

function plotLines(points, color, name, LineWidth)
for i = 1:(size(points, 2)/2)
    line(points(1, 2*i-1:2*i), points(2, 2*i-1:2*i), 'Color', color, 'LineWidth', LineWidth);
    midPoint = mean(points(:, 2*i-1:2*i), 2);
    text(midPoint(1), midPoint(2), [name num2str(i)], 'Color', color, 'FontSize', 11);
end
end

function plotConicParams(center, axes, angle, color)
%% Plot the axes
roationMatrix = [cos(angle), -sin(angle); sin(angle), cos(angle)];
axesVectors = roationMatrix * diag(axes);
quiver(center(1), center(2), axesVectors(1, 1), axesVectors(2, 1), 'Color', color, 'LineWidth', 2);
quiver(center(1), center(2), axesVectors(1, 2), axesVectors(2, 2), 'Color', color, 'LineWidth', 2);

%% Plot the center
plot(center(1), center(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end