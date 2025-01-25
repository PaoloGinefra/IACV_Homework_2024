close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Load the lines for visualization purposes
load('../lines.mat');

%% Load the calibration matrix
load('../K_metric.mat');

new_r = cross(rs(:, 2), rs(:, 1));
new_r = new_r / norm(new_r) * norm(rs(:, 1));

new_rs = [rs(:, 1), new_r, rs(:, 3)];

%% Plot hte original image
figure;
imshow(im);

%% Plot in a 3d scatter plot rs(:, 3)
figure;
plot3(rs(1, 3), rs(2, 3), rs(3, 3), 'xb');
hold on;
% plot the other columns as arrows from the 3rd
quiver3(rs(1, 3), rs(2, 3), rs(3, 3), rs(1, 1), rs(2, 1), rs(3, 1), 'Color', 'r');
quiver3(rs(1, 3), rs(2, 3), rs(3, 3), rs(1, 2), rs(2, 2), rs(3, 2), 'Color', 'g');
quiver3(rs(1, 3), rs(2, 3), rs(3, 3), new_r(1), new_r(2), new_r(3), 'Color', 'b');




H = K * new_rs;

H_rect = inv(H);

%% warp the image
tform = projective2d(H_rect');
im_warped = imwarp(im, tform);

% plot the image
figure;
imshow(im_warped);
title('Rectified Image');
impixelinfo;


%% Try for s

rs_for_S = [rs(:, 1), rs(:, 2), rs(:, 3) + new_r / 2];
H = K * rs_for_S;

H_rect = inv(H);

%% warp the image
tform = projective2d(H_rect');
im_warped = imwarp(im, tform);

% plot the image
figure;
imshow(im_warped);
title('Rectified Image for s');