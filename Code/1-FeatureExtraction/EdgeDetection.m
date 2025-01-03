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

%% Apply Cannys edge detection
threshold = [0.05, 0.15];
sigma = 5;
im_edges = edge(im_gray, 'canny', threshold, sigma);

% mask the edges
im_edges = im_edges .* image_mask;

%% Display the results
figure;
imshow(im_edges);
title(['Canny Edges with mask - Threshold: ', num2str(threshold), ' Sigma: ', num2str(sigma)]);