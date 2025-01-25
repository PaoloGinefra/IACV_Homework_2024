close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Load the lines for visualization purposes
load('../2.0-ManualLineExtraction/lines.mat');

%% Load the metric rectification Homography and vars
load('../2.2-Metric_Rectification/H_metric.mat');

%% Load the calibration matrix
load('../2.3-Intrinsic_Calibration/K.mat');

%% Height Calculation method 1
height = (1 - l1_length/l2_length) * rs(:, 3)' * r3;
disp('Height');
disp(height);

save('height.mat', 'height');

new_r = cross(rs(:, 2), rs(:, 1));
new_r = new_r / norm(new_r) * norm(rs(:, 1));

new_rs = [rs(:, 1), new_r, rs(:, 3)];

H = K * new_rs;

H_rect = [10, 0, 0; 0 10 0; 0, 0, 1] * inv(H);

%% warp the image
tform = projective2d(H_rect');
im_warped = imwarp(im, tform);

% plot the image
figure;
imshow(im_warped);
title('Rectified Image');