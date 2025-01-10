close all;
clc;
clear;

%% Load Image
im = imread('../../Assignment/Homework Image.jpg');

crop_top = 220;
crop_bottom = 823;
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

%% Plot histogram per feature
% figure;
% for i = 1:size(features, 1)
%     subplot(2, 3, i);
%     histogram(features(i, :), 100);
%     title(['Feature ', num2str(i)]);
% end


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
imagesc(reshape(features_pca(1, :), size(im, 1), size(im, 2)));
title('First Principal Component');
impixelinfo;
%set gray colormap
colormap(gray);

% Plot the image in the second principal component
figure;
imagesc(reshape(features_pca(2, :), size(im, 1), size(im, 2)));
title('Second Principal Component');
impixelinfo;
%set gray colormap
colormap(gray);

% Plot the image in the third principal component
figure;
imagesc(reshape(features_pca(3, :), size(im, 1), size(im, 2)));
title('Third Principal Component');
impixelinfo;
%set gray colormap
colormap(gray);

% Mask
threshold_upper = 0;
threshold_lower = -0.6;
mask = features_pca(1, :) < threshold_upper & features_pca(1, :) > threshold_lower;


% Plot the mask
figure;
imagesc(reshape(mask, size(im, 1), size(im, 2)));
title('Mask');

% Convolution filtering
%gaussian filter
sigma = 15;
kernel_size = 2 * sigma;
gaussian_filter = fspecial('gaussian', kernel_size, sigma);
filter_max = max(gaussian_filter(:));

gaussian_filter = 2  * gaussian_filter - filter_max;

% Apply the filter
source = reshape(features_pca(2, :), size(im, 1), size(im, 2));
filtered_image = conv2(source, gaussian_filter, 'same');

% Plot the filtered image
figure;
imagesc(filtered_image);
title('Filtered Image');
impixelinfo;

% Threshold the filtered image
threshold = -0.05;
filtered_mask = filtered_image < threshold;

% Plot the filtered mask
figure;
imagesc(filtered_mask);
title('Filtered Mask');





