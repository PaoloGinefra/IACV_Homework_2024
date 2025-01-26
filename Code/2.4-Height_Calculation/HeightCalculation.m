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