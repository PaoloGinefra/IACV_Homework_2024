close all;
clc;
clear;

%% Load Image
im = imread('../../Assignment/Homework Image.jpg');

crop_top = 220;
crop_bottom = 800;
crop_left = 285;
crop_right = 1462;
im = im(crop_top+1:crop_bottom, crop_left+1:crop_right, :);

figure;
imshow(im);
impixelinfo;

%% Convert to hsv
im_hsv = rgb2hsv(im);

%% reshape the image to be 3, cols*rows
RGB_features = reshape(double(im) / 255, size(im, 1)*size(im, 2), 3)';
HSV_features = reshape(im_hsv, size(im_hsv, 1)*size(im_hsv, 2), 3)';

% features = [RGB_features(1:2, :); HSV_features(2:2, :)];
features = [RGB_features; HSV_features(2:3, :)];


%% Perform PCA

% Subtract the mean
mean_features = mean(features, 2);
features_centered = features - mean_features;


% Compute the covariance matrix
covariance_matrix = cov(features_centered');

% Plot the covariance matrix
figure;
imagesc(covariance_matrix);
title('Covariance Matrix');
colorbar();

% Perform SVD
[U, S, V] = svd(features_centered, "econ");

% Plot the singular values
figure;
subplot(1, 3, 1);
semilogy(diag(S));

% Plot the cumulative sum of singular values
subplot(1, 3, 2);
plot(cumsum(diag(S))/sum(diag(S)));

% Plot the explained variance
subplot(1, 3, 3);
plot(cumsum(power(diag(S), 2))/sum(power(diag(S), 2)));

% Project the features onto the principal components
num_components = 3;
features_pca = U(:, 1:num_components)' * features_centered;

% Plot the PCA features
figure;
scatter(features_pca(1, :), features_pca(2, :), 1, 'filled');
title('PCA Features');
xlabel('Principal Component 1');
ylabel('Principal Component 2');

% Plot the PCA features in 3d, coloring using the original features
figure;
scatter3(features_pca(1, :), features_pca(2, :), features_pca(3, :), 2, RGB_features', 'filled');
title('PCA Features');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
zlabel('Principal Component 3');

% Plot the image in the first principal component
figure;
subplot(3, 1, 1);
imagesc(reshape(features_pca(1, :), size(im, 1), size(im, 2)));
title('First Principal Component');
axis equal;
impixelinfo;
%set gray colormap
colormap(gray);

% Plot the image in the second principal component
subplot(3, 1, 2);
imagesc(reshape(features_pca(2, :), size(im, 1), size(im, 2)));
title('Second Principal Component');
axis equal;
impixelinfo;
%set gray colormap
colormap(gray);

% Plot the image in the third principal component
subplot(3, 1, 3);
imagesc(reshape(features_pca(3, :), size(im, 1), size(im, 2)));
title('Third Principal Component');
axis equal;
impixelinfo;
%set gray colormap
colormap(gray);

% Canny edge detection on first principal component
threshold = [0.0, 0.08];
sigma = 1.5;
im_edges = edge(reshape(features_pca(1, :), size(im, 1), size(im, 2)), 'canny', threshold, sigma);

% Plot the edges
figure;
imshow(im_edges);
title('Canny Edges on First Principal Component');

% blur the edges
sigma = 20;
im_edges_blurred = imgaussfilt(double(im_edges), sigma);

figure;
subplot(2, 1, 1);
imagesc(im_edges_blurred);
axis equal;
title('Blurred Edges');

% threshold the edges
threshold = 0.065;
im_edges_mask_pure = im_edges_blurred < threshold;

subplot(2, 1, 2);
imagesc(im_edges_mask_pure);
axis equal;
title('Thresholded Edges');

% mask the edges
im_edges = im_edges_mask_pure .* im_edges;

figure;
imshow(im_edges);
axis equal;
title('Masked Edges');

% apply hough transform
[H, theta, rho] = hough(im_edges, "Theta", -90:1:89.9, "RhoResolution", 6);
peaks = houghpeaks(H, 50, "Threshold", 0.3 * max(H(:)), "NHoodSize", [11, 11]);
lines = houghlines(im_edges, theta, rho, peaks, "FillGap", 50, "MinLength", 150);

% Plot the lines
figure;
imshow(im_edges);
hold on;
for i = 1:length(lines)
    xy = [lines(i).point1; lines(i).point2];
    plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'red');
end
title('Hough Transform Lines');


