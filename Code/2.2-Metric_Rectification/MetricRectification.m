close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Load the lines
% Loads: l_points, l_lines, m_points, m_lines, h_points, h_lines
load('../2.0-ManualLineExtraction/lines.mat');

%% load the line at infinity
% Loads: l_h_inf_prime
load('../2.1-Vanishing_Line_Extraction/VanishingLineHorizontalPlane.mat');

%% Load the extracted conic
load('../1-FeatureExtraction/ConicExtraction/ExtractedConic.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Affine Rectification                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the homography for the affine rectification
l_h_inf_prime = l_h_inf_prime / l_h_inf_prime(3);

% Added a scaling factor of 0.1 to make the rectified image more managable
H_aff = [.1, 0, 0; 0, .1, 0; l_h_inf_prime'];
H_aff_inv = inv(H_aff);

%% Sanity check the homography
assert(all(abs(warpLine(l_h_inf_prime, H_aff_inv)' - [0, 0, 1]) < 1e-6));
disp('Sanity Check [Affine Rectification] - Line at infinity after homography is at infinity: ✅');

%% Apply the homography to the image
im_warped = warpImage(im, H_aff);

%% Compute the image of the lines
m_lines_aff = warpLine(m_lines, H_aff_inv);
m_lines_aff = m_lines_aff ./ m_lines_aff(3, :);
m_points_aff = warpPoint(m_points, H_aff);

l_lines_aff = warpLine(l_lines, H_aff_inv);
l_lines_aff = l_lines_aff ./ l_lines_aff(3, :);
l_points_aff = warpPoint(l_points, H_aff);

%% Sanity check: parallel lines
m_lines_aff_angles = rad2deg(atan2(m_lines_aff(2, :), m_lines_aff(1, :)));
l_lines_aff_angles = rad2deg(atan2(l_lines_aff(2, :), l_lines_aff(1, :)));
assert(all(abs(m_lines_aff_angles - [m_lines_aff_angles(2:end), m_lines_aff_angles(1)]) < 1));
assert(all(abs(l_lines_aff_angles - [l_lines_aff_angles(2:end), l_lines_aff_angles(1)]) < 1));
disp('Sanity Check [Affine Rectification] - Parallel lines after affine rectification have the same direction up to 1 deg: ✅');

%% Compute conic Errors
conicErrors = computeConicErrors(C, im);
[center, axes, angle] = paramsFromHomogenousConic_separated(C);

%% Apply the homography
C_aff = warpConic(C, H_aff_inv);
C_aff_par = paramsFromHomogenousConic(C_aff);

%% Compute the conic errors after the affine rectification
C_aff_errors = computeConicErrors(C_aff, im_warped);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Metric rectification                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% convert the conic coefficient to geometric parameters
[center_aff, axes_aff, angle_aff] = paramsFromHomogenousConic_separated(C_aff);

% Rotation
U = [cos(angle_aff), -sin(angle_aff); sin(angle_aff), cos(angle_aff);];

% Rescaling the axis to make them equal
a = axes_aff(1);
b = axes_aff(2);
S = diag([1, a/b]);
K = U*S*U';
H_2 = [K, zeros(2, 1); 0, 0, 1];

% Compute the metric rectification homography
H_metric = H_2 * H_aff;
H_metric_inv = inv(H_metric);

%% Apply rotation to make l lines paralle to the x-axis
% Compute the warped lines
l_lines_metric = warpLine(l_lines, H_metric_inv);

% Compute the average angle of the lines
angles = rad2deg(atan2(l_lines_metric(2, :), l_lines_metric(1, :)));
avg_rotation = -mean(angles) + 90;

% Compute the rotation matrix
R = [cosd(avg_rotation), -sind(avg_rotation); sind(avg_rotation), cosd(avg_rotation)];
H_rotation = [R, zeros(2, 1); 0, 0, 1];

% Apply the rotation to the homography
H_metric = H_rotation * H_metric;
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
assert(all(norm(m_lines_metric(1:2, :)' * l_lines_metric(1:2, :)) < 1e-3));
disp('Sanity Check [Metric Rectification] - The rectified lines are orthogonal: ✅');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  Depth                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute depth m using the length of l1 and the average length of the m lines
m_lines_length = vecnorm(m_points_metric(:, 1:2:end) - m_points_metric(:, 2:2:end));

l1_length = vecnorm(l_points_metric(:, 1) - l_points_metric(:, 2));
l2_length = vecnorm(l_points_metric(:, 3) - l_points_metric(:, 4));

average_m_length = mean(m_lines_length);
depth_m = average_m_length / l1_length;

disp(['Depth of m: ' num2str(depth_m)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 Plotting                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%% Save the matric rectification homography
% H_metric = H_metric/l1_legnth;
save('./H_metric.mat', 'H_metric', 'l1_length', "l2_length", 'l_points_metric', 'm_points_metric', 'depth_m');


%% Pick 12 points on the curve S
PICK_POINTS = false;
n_points = 12;

if PICK_POINTS
    figure;
    imshow(im_metric);
    title(['Pick ' num2str(n_points) ' points on the curve S']);
    hold on;
    [x, y] = ginput(n_points);
    S_points = [x, y];
    plot(S_points(:, 1), S_points(:, 2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    % label points
    for i = 1:n_points
        text(S_points(i, 1), S_points(i, 2), ['S' num2str(i)], 'Color', 'r', 'FontSize', 11);
    end
    hold off;
    save('./S_points.mat', 'S_points');
end


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

h_warped = ceil(max(corners_warped(2, :)) - min(corners_warped(2, :)));
w_warped = ceil(max(corners_warped(1, :)) - min(corners_warped(1, :)));

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