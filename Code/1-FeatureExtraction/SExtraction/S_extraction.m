close all;
clc;
clear;

%% Load the image
im = imread('../../../Assignment/Homework Image.jpg');

figure;
imshow(im);
title('Original Image');
impixelinfo;

%% Cut out the region of interest
% Define the region of interest
BL_corner = [830, 546];
UR_corner = [952, 510];

% Cut out the region of interest
im = im(UR_corner(2):BL_corner(2), BL_corner(1):UR_corner(1), :);
figure;
subplot(3, 1, 1);
imshow(im);
title('Cropped Image');

%% Canny Edge Detection
% Convert image to grayscale
im_gray = rgb2gray(im);

subplot(3, 1, 2);
imshow(im_gray);
title('Grayscale Image');

% Apply Cannys edge detection
threshold = [0.2, 0.8];
sigma = 1;
im_edges = edge(im_gray, 'canny', threshold, sigma);

% Display the results
subplot(3, 1, 3);
imshow(im_edges);
title(['Canny Edges - Threshold: ', num2str(threshold), ' Sigma: ', num2str(sigma)]);

% Save the plot
saveas(gcf, 'CannyEdges.png');