close all;
clc;
clear;

%% Load the image
im = imread('../Assignment/Homework Image.jpg');

%% load the line at infinity
load('./VanishingLineHorizontalPlane.mat');

%% load conic
load('./CircleC.mat');

%% Compute the homography
H = [-0.1, 0, 0; 0, -0.1, 0; l_h_inf_prime'];

%% Apply the homography
im_warped = warpImage(im, H);

%% Sanity check the homography
disp('Sanity Check: line at infinity after homography');
disp(inv(H)' * l_h_inf_prime);

%% Show the images
figure;
subplot(1, 2, 1);
imshow(im);
title('Original Image');

% Plot The conic
conicPlot = zeros(size(im, 1), size(im, 2));

for i = 1:size(im, 1)
    for j = 1:size(im, 2)
        point = [j, i, 1];
        conicPlot(i, j) = point * C * point';
    end
end

%% Apply the homography
H_inv = inv(H);
C_warped = H_inv' * C * H_inv;
C_warped = C_warped ./ C_warped(3, 3);

disp(C_warped);

disp('Conic Warped');

%% Plot the conic
conicPlotWarped = zeros(size(im_warped, 1), size(im_warped, 2));

for i = 1:size(im_warped, 1)
    for j = 1:size(im_warped, 2)
        point = [j, i, 1];
        conicPlotWarped(i, j) = point * C_warped * point';
    end
end

%% convert the conic coefficient to geometric parameters
par_geo = AtoG([C_warped(1,1),2*C_warped(1,2),C_warped(2,2),2*C_warped(1,3),2*C_warped(2,3),C_warped(3,3)]);
center = par_geo(1:2);
axes = par_geo(3:4);
angle = par_geo(5);


disp('Conic Plot Computed');

subplot(1, 2, 1);
imshow(im);
hold on;
contour(conicPlot, [0, 0], 'r', 'LineWidth', 2);
% contour(conicPlotWarped, [0, 0], 'r', 'LineWidth', 2);
title('Original image + Conic');
hold off;

subplot(1, 2, 2);
imshow(im_warped);
hold on;
contour(conicPlotWarped, [0, 0], 'r', 'LineWidth', 2);
plot(center(1),center(2),'ro','Markersize',20);
title('Image + affine rectification');
hold off;


%% Metric rectification
% Rotation
U = [cos(angle), -sin(angle); sin(angle), cos(angle);];

% rescaling the axis to make them equal
a = axes(1);
b = axes(2);
S = diag([1, a/b]);
K = U*S*U';
H_2 = [K, zeros(2, 1); 0, 0, 1];

H_metric = H_2 * H;

%% Apply the metric rectification
im_metric = warpImage(im, H_metric);

%% Compute the image of the conic
H_metric_inv = inv(H_metric);
C_metric = H_metric_inv' * C * H_metric_inv;

disp('Conic Metric');
disp(C_metric);

%% Plot the conic
conicPlotMetric = zeros(size(im, 1), size(im, 2));

for i = 1:size(im, 1)
    for j = 1:size(im, 2)
        point = [j, i, 1];
        conicPlotMetric(i, j) = point * C_metric * point';
    end
end

disp('Conic Plot Metric Computed');

%% Compute the images of the lines
load("./lines.mat")
m_lines_metric = H_metric_inv' * m_lines;
m_lines_metric = m_lines_metric ./ vecnorm(m_lines_metric);
m_points_metric = H_metric * [m_points; ones(1, size(m_points, 2))];
m_points_metric_euclidian = m_points_metric(1:2, :) ./ m_points_metric(3, :);

l_lines_metric = H_metric_inv' * l_lines;
l_lines_metric = l_lines_metric ./ vecnorm(l_lines_metric);
l_points_metric = H_metric * [l_points; ones(1, size(l_points, 2))];
l_points_metric_euclidian = l_points_metric(1:2, :) ./ l_points_metric(3, :);

%% Sanity check: othogonality of rectified lines
disp('Sanity Check: Orthogonality of rectified lines');
disp(m_lines_metric' * l_lines_metric);

%% Show the images
figure;
subplot(1, 2, 1);
imshow(im);
title('Original Image');

subplot(1, 2, 2);
imshow(im_metric);
hold on;
contour(conicPlotMetric, [0, 0], 'r', 'LineWidth', 2);

for i = 1:size(m_lines_metric, 2)
    line(m_points_metric_euclidian(1, 2*i-1:2*i), m_points_metric_euclidian(2, 2*i-1:2*i), 'Color', 'g', 'LineWidth', 2);
end

for i = 1:size(l_lines_metric, 2)
    line(l_points_metric_euclidian(1, 2*i-1:2*i), l_points_metric_euclidian(2, 2*i-1:2*i), 'Color', 'r', 'LineWidth', 2);
end
hold off
title('Image + metric rectification');


function [imageWarped] = warpImage(image, H)
h = size(image, 1);
w = size(image, 2);

corners =  [
    0, w, 0, w;
    0, 0, h, h;
    ones(1, 4)
    ];

corners_warped = H * corners;
corners_warped = corners_warped ./ corners_warped(3, :);

disp('corners_warped');
disp(corners_warped);

w_warped = ceil(max(corners_warped(1, :)) - min(corners_warped(1, :)));
h_warped = ceil(max(corners_warped(2, :)) - min(corners_warped(2, :)));
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