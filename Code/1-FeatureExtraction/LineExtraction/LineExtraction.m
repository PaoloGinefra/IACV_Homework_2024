close all;
clc;
clear;

%% Load Image
im = imread('../../../Assignment/Homework Image.jpg');

%% Load mask
load('../../1-FeatureExtraction/image_mask.mat');

%% Convert to hsv
im_hsv = rgb2hsv(im);

%% Convert to grayscale
im_gray = im(:,:,2); %rgb2gray(im);

%% Plot the grayscale image
figure;
imshow(im_gray);
title('Grayscale Image');

%% Canny edge detection
threshold = [0, 0.04];
sigma = 1;
im_edges = edge(im_gray, 'canny', threshold, sigma);

%mask the edges
% im_edges = im_edges .* image_mask;

%% Display the results
figure;
imshow(im_edges);
title(['Canny Edge Detection with threshold = [', num2str(threshold(1)), ', ', num2str(threshold(2)), '] and sigma = ', num2str(sigma)]);

% Blur the edges
sigma = 8;
im_edges_blured = imgaussfilt(double(im_edges), sigma);

% threshold the edges
threshold = 0.11;
im_edges_mask = im_edges_blured < threshold;

figure;
subplot(1, 2, 1);
imshow(im_edges_blured);
title(['Blurred Edges with sigma = ', num2str(sigma),]);

subplot(1, 2, 2);
imshow(im_edges_mask);
title(['Thresholded Edges with threshold = ', num2str(threshold)]);
impixelinfo;

% Apply the mask
im_edges_mask = im_edges_mask .* image_mask;
im_edges = im_edges .* im_edges_mask;

figure;
imshow(im_edges);
title('Canny Edge Detection');

%% Compute the gradients of the edges
sigma = 2;
im_gray_smoothed = imgaussfilt(im_gray, sigma);
[Gx, Gy] = imgradientxy(im_gray_smoothed);

%% Compute the gradient direction
[Gmag, Gdir] = imgradient(Gx, Gy);

% plot the gradient direction as hue
Gdir_color_hsv = zeros(size(im_edges, 1), size(im_edges, 2), 3);
Gdir_color_hsv(:, :, 1) = (Gdir + 180) / 360;
Gdir_color_hsv(:, :, 2) = ones(size(im_edges, 1), size(im_edges, 2));
Gdir_color_hsv(:, :, 3) = ones(size(im_edges, 1), size(im_edges, 2));

Gdir_color = hsv2rgb(Gdir_color_hsv);

Gdir_color = Gdir_color .* im_edges;

figure;
imshow(Gdir_color);
title('Gradient Direction');

%% Find the lines
[y, x] = find(im_edges> 0.7);

points = [x, y, ones(size(x, 1), 1)];


%plot the points
figure;
imshow(im_edges);
hold on;
plot(points(:, 1), points(:, 2), 'r.', 'MarkerSize', 10);
title('Edge Points');


normals = zeros(size(points, 1), 2);
for i = 1:size(points, 1)
    g = Gdir(points(i, 2), points(i, 1));
    normals(i, 1) = cosd(g);
    normals(i, 2) = sind(g);
end


pointsAtInfinity = [normals(:, 2), normals(:, 1), zeros(size(normals, 1), 1)];

% normalize the points with image dimensions
intercept_weight = 1;
points_normalized = points ./ [size(im, 2)/intercept_weight, size(im, 1)/intercept_weight, 1];

lines = cross(points_normalized, pointsAtInfinity);

%flip lines with negative third element
lines = lines .* sign(lines(:, 3));

%normalize the lines
% lines = lines ./ vecnorm(lines, 2, 2);

% plot the lines as points in a 3d scatter plot
figure;
hold on;
scatter3(lines(:, 1), lines(:, 2), lines(:, 3), 'r.');

%% Use DBSCAN to cluster the lines
epsilon = 0.02;
minPts = 40;
[labels, numClusters] = dbscan(lines, epsilon, minPts);

n_clusters = max(labels);
disp(['Number of clusters: ', num2str(n_clusters)]);

colors = rand(n_clusters, 3);

% plot the clusters
figure;
hold on;
for i = 1:n_clusters
    cluster = lines(labels == i, :);
    scatter3(cluster(:, 1), cluster(:, 2), cluster(:, 3), 'filled', 'MarkerFaceColor', colors(i, :));
end
title(['DBSCAN Clustering with epsilon = ', num2str(epsilon), ' and minPts = ', num2str(minPts)]);
grid on;
hold off;


% compute centroids of the clusters
centroids = zeros(n_clusters, 3);
for i = 1:n_clusters
    cluster = lines(labels == i, :);
    centroids(i, :) = mean(cluster, 1);
end

% plot the centroids as lines in the image
figure;
imshow(im);
hold on;
for i = 1:n_clusters
    line = centroids(i, :);
    line = line ./ line(3);
    plot([1, size(im, 2)], [(-line(3) / line(2))*size(im, 1)/intercept_weight, (-line(3) - line(1)) / line(2) * size(im, 1)/intercept_weight], 'Color', colors(i, :), 'LineWidth', 2);
end


%% Use the Hough Transform to find the lines
[H, theta, rho] = hough(im_edges, "Theta", -90:0.01:89.9);

peaks = houghpeaks(H, 50, 'Threshold', 0.2 * max(H(:)));

lines = houghlines(im_edges, theta, rho, peaks, 'FillGap', 500, 'MinLength', 100);

% plot the lines
figure;
imshow(im);
hold on;
for i = 1:length(lines)
    xy = [lines(i).point1; lines(i).point2];
    plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'green');
end
title('Hough Transform Lines');
