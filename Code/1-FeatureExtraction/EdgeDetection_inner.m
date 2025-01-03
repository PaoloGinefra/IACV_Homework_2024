close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Load the mask
load('./image_mask.mat');

im = double(im) .* image_mask;
im = uint8(im);

%% Convert image to hsv
im_hsv = rgb2hsv(im);

%% Convert the image to grayscale
im_gray = im(:, :, 1);%%rgb2gray(im);

%% Rescale the image
min = 60;
max = 120;
im_gray = imadjust(im_gray, [min/255, max/255], [0, 1]);

figure;
imshow(im_gray);
title('Grayscale Image');


%% Apply Cannys edge detection
threshold = [0.05, 0.1];
sigma = 5;
im_edges = edge(im_gray, 'canny', threshold, sigma);

% mask the edges
im_edges = im_edges .* image_mask;

%% Display the results
figure;
imshow(im_edges);
title(['Canny Edges with mask - Threshold: ', num2str(threshold), ' Sigma: ', num2str(sigma)]);

%% Overlay edges on the original image
im_edges_rgb = repmat(im_edges, [1, 1, 3]);
im_edges_rgb(:, :, 2:3) = 0;

im_edges_rgb = uint8(im_edges_rgb) .* 255;

im_overlay = im;
im_overlay(im_edges_rgb == 255) = 255;

figure;
imshow(im_overlay);
title('Overlay of edges on the original image');