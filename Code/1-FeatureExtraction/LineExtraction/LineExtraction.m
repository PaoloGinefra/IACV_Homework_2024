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
sigma = 2;
im_edges = edge(im_gray, 'canny', threshold, sigma);
im_edges = im_edges .* image_mask;

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
im_edges_mask_pure = im_edges_blured < threshold;

% Erode the mask
erosion_size = 10;
sd = strel('disk', erosion_size);
im_edges_mask = imdilate(im_edges_mask_pure, sd);
se = strel('disk', erosion_size + 15);
im_edges_mask = imerode(im_edges_mask, se);

figure;
subplot(1, 3, 1);
imshow(im_edges_blured);
title(['Blurred Edges with sigma = ', num2str(sigma),]);

subplot(1, 3, 2);
imshow(im_edges_mask_pure);
title(['Thresholded Edges with threshold = ', num2str(threshold)]);

subplot(1, 3, 3);
imshow(im_edges_mask);
title(['Eroded Edges with erosion size = ', num2str(erosion_size)]);

impixelinfo;

% Apply the mask
im_edges_mask = im_edges_mask .* image_mask;
im_edges = im_edges .* im_edges_mask;

im_edges = im_edges > 0.5;

figure;
imshow(im_edges);
title('Canny Edge Detection');

%% Use the Hough Transform to find the lines
[H, theta, rho] = hough(im_edges, "Theta", -90:1:89.9);

peaks = houghpeaks(H, 50, 'Threshold', 0.3 * max(H(:)), 'NHoodSize', [21, 21]);

lines = houghlines(im_edges, theta, rho, peaks, 'FillGap', 150, 'MinLength', 100);

% plot the lines
figure;
imshow(im);
hold on;
for i = 1:length(lines)
    xy = [lines(i).point1; lines(i).point2];
    plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'green');
end
title('Hough Transform Lines');
