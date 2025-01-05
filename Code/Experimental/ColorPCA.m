close all;
clc;
clear;

%% Load Image
im = imread('../../Assignment/Homework Image.jpg');

figure;
imshow(im);

%% Convert to hsv
im_hsv = rgb2hsv(im);

%% reshape the image to be 3, cols*rows
RGB_features = reshape(double(im) / 255, size(im, 1)*size(im, 2), 3)';
HSV_features = reshape(im_hsv, size(im_hsv, 1)*size(im_hsv, 2), 3)';

features = [RGB_features(1:2, :); HSV_features(2:2, :)];

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

% Plot a 2d histogram of the PCA features
figure;




