close all;
clc;
clear;

%% Load Image
im = imread('../../Assignment/Homework Image.jpg');

%% Load mask
load('../1-FeatureExtraction/image_mask.mat')

%% Convert to hsv
im_hsv = rgb2hsv(im);

%% Convert to grayscale
im_gray = im_hsv(:, :, 2); %rgb2gray(im);

%% Smooth the image
sigma = 1;
im_gray = imgaussfilt(im_gray, sigma);

%% Compute the x and y gradients
[Gx, Gy] = imgradientxy(im_gray);

%% Compute the gradient magnitude and direction
[Gmag, Gdir] = imgradient(Gx, Gy);

% Smooth gdir
% sigma_gdir = 1;
% Gdir = imgaussfilt(Gdir, sigma_gdir);

Gdir_color_hsv = zeros(size(im_gray, 1), size(im_gray, 2), 3);
Gdir_color_hsv(:, :, 1) = (Gdir + 180) / 360;
Gdir_color_hsv(:, :, 2) = ones(size(im_gray, 1), size(im_gray, 2));
Gdir_color_hsv(:, :, 3) = ones(size(im_gray, 1), size(im_gray, 2));

Gdir_color = hsv2rgb(Gdir_color_hsv);

figure;
imshow(Gdir_color);
title('Gradient Direction');

%% Display the results
figure;
subplot(2, 2, 1);
imshow(im_gray);
title('Original Image');

subplot(2, 2, 2);
imshow(Gmag, [0, 1]);
title('Gradient Magnitude');

subplot(2, 2, 3);
imshow(Gdir, []);
colormap jet;
colorbar;
title('Gradient Direction');

subplot(2, 2, 4);
imshowpair(Gx, Gy, 'montage');
title('Gradient X and Y');

impixelinfo;

%% Canny edge detection
threshold = [0.05, 0.2];
sigma = 1.5;
im_edges = edge(im_gray, 'canny', threshold, sigma);

% mask the edges
im_edges = im_edges .* image_mask;

figure;
imshow(im_edges);
title(['Canny Edges with mask - Threshold: ', num2str(threshold), ' Sigma: ', num2str(sigma)]);

% smooth the edges
sigma_edges = 8;
im_edges_smoothed = imgaussfilt(im_edges, sigma_edges);

figure;
imshow(im_edges_smoothed);
title(['Smoothed Edges - Sigma: ', num2str(sigma_edges)]);
impixelinfo;

im_edges_smoothed_threshold = 0.1;
im_edges_smoothed(im_edges_smoothed >= im_edges_smoothed_threshold) = 1;
im_edges_smoothed(im_edges_smoothed < im_edges_smoothed_threshold) = 0;

%dilate the mask
se = strel('disk', 5);
im_edges_smoothed = imdilate(im_edges_smoothed, se);

figure;
imshow(im_edges_smoothed);
title(['Thresholded Smoothed Edges - Threshold: ', num2str(im_edges_smoothed_threshold)]);

im_edges = im_edges .* (1 - im_edges_smoothed);

figure;
imshow(im_edges);
title('Edges with smoothed edges removed');
impixelinfo;

% Clamp the gradient magnitude to [0, 1]
% mask the gradient magnitude using the image mask
Gmag_clamped = Gmag .* image_mask;
Gmag_threshold = 0.4;
Gmag_clamped(Gmag_clamped >= Gmag_threshold) = 1;
Gmag_clamped(Gmag_clamped < Gmag_threshold) = 0;

figure;
imshow(Gmag_clamped);
title('Clamped Gradient Magnitude');

% mask the gradient dir with the gradient magnitude
% Gdir_mask = Gmag_clamped;
Gdir_mask = im_edges;
Gdir_masked = Gdir_color;
Gdir_masked(:, :, 1) = Gdir_masked(:, :, 1) .* Gdir_mask;
Gdir_masked(:, :, 2) = Gdir_masked(:, :, 2) .* Gdir_mask;
Gdir_masked(:, :, 3) = Gdir_masked(:, :, 3) .* Gdir_mask;

figure;
imshow(Gdir_masked);
title('Masked Gradient Direction');

% Select the gradient direction based on the gradient magnitude
Gdir_rad = Gdir * pi / 180;
Gdir_x = cos(Gdir_rad);
Gdir_y = -sin(Gdir_rad);

[ys, xs] = find(im_edges > 0.5);
points = [xs / size(im, 2), ys/size(im, 1), ones(size(xs))]';

Gdir_x_selected = Gdir_x(sub2ind(size(Gdir_x), ys, xs));
Gdir_y_selected = Gdir_y(sub2ind(size(Gdir_y), ys, xs));

pointsAtInfinity = [-Gdir_y_selected, Gdir_x_selected, zeros(size(Gdir_x_selected))]';

lines = cross(points, pointsAtInfinity);
%lines = lines ./ size(im, 2);
lines = lines .* sign(lines(3, :));
lines = lines ./ sqrt(lines(1, :).^2 + lines(2, :).^2 + lines(3, :).^2);

% Plot the lines
figure;
imshow(im);
hold on;
randomIndices = randperm(size(lines, 2));
for j = 1:10
    i = randomIndices(j);
    plot([1, size(im, 2)], [(-lines(3, i) / lines(2, i))*size(im, 1), (-lines(3, i) - lines(1, i)) / lines(2, i) * size(im, 1)], 'Color', 'r');
    % Use Gdir_color as the color
    plot(xs(i), ys(i), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', Gdir_color(ys(i), xs(i), :));
    % draw an arrow in the direction of the gradient
    quiver(xs(i), ys(i), Gdir_x_selected(i) * 10, Gdir_y_selected(i) * 10, 'Color', 'r');
end
hold off;
title('Selected Gradient Direction');

%% Make a 3d scatterplot of the lines array
figure;
hold on;
% use random indices
randomLines = randomIndices(1:100);
% Use Gdir_color as the color
for k = 1:length(randomLines)
    scatter3(lines(1, randomLines(k)), lines(2, randomLines(k)), lines(3, randomLines(k)), 'MarkerEdgeColor', squeeze(Gdir_color(ys(randomLines(k)), xs(randomLines(k)), :))');
end
xlabel('a');
ylabel('b');
zlabel('c');

% add grid
grid on;

hold off;



%% Clustering of the lines
% Use DbScan to cluster the lines
epsilon = 0.008;
% density = 1000000;
minPts = 40;
disp(['Epsilon: ', num2str(epsilon)]);
disp(['MinPts: ', num2str(minPts)]);
[clusterIdx, isNoise] = dbscan(lines', epsilon, minPts);

disp(['Number of clusters: ', num2str(max(clusterIdx))]);

% Compute the centroid of each cluster
centroids = zeros(3, max(clusterIdx));
clusterColors = zeros(max(clusterIdx), 3);
for i = 1:max(clusterIdx)
    cluster = lines(:, clusterIdx == i);
    centroids(:, i) = mean(cluster, 2);
    clusterColors(i, :) = rand(1, 3);
end

% Plot the clusters
figure();
hold on;
for i = 1:max(clusterIdx)
    cluster = lines(:, clusterIdx == i);
    scatter3(cluster(1, :), cluster(2, :), cluster(3, :), 'MarkerEdgeColor', clusterColors(i, :));
end
hold off;
grid on;

% Plot the centroids as lines on the image
figure();
imshow(im);
hold on;
for i = 1:max(clusterIdx)
    plot([1, size(im, 2)], [(-centroids(3, i) / centroids(2, i))*size(im, 1), (-centroids(3, i) - centroids(1, i)) / centroids(2, i) * size(im, 1)], 'Color', clusterColors(i, :));
end
hold off;

