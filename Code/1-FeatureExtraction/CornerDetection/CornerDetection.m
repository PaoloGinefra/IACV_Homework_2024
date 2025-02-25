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
im_gray = im_gray .* uint8(image_mask);

% HARRIS ALGORITHM: find corner features in an image

% Derivative masks
dx = [-1 0 1;
    -1 0 1;
    -1 0 1];
dy = dx';

% Image derivatives : applying the filter through convolution
Ix = conv2(im_gray, dx, 'same');
Iy = conv2(im_gray, dy, 'same');

% Define sigma
sigma = 3; % You can adjust this value as needed

% Gaussian filter
g = fspecial('gaussian',max(1,fix(3*sigma)+1), sigma);

% applying the gaussian filter to images
Ix2 = conv2(Ix.^2, g, 'same');
Iy2 = conv2(Iy.^2, g, 'same');
Ixy = conv2(Ix.*Iy, g, 'same');

% combining filtered images
cm = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps);

% Define the margin for the border
border_margin=5;

% Set to 0 near the boundaries
cm(1:border_margin,:)=0;
cm(end-border_margin:end,:)=0;
cm(:,end-border_margin:end)=0;
cm(:,1:border_margin)=0;

% Threshold the cim
T=mean(cm(:));
CIM=cm;
CIM(cm<T)=0;

% perform nonlocal maximum suppression on the thresholded measure
support=true(11);

% compute maximum over a square neighbor of size 11 x 11
maxima=ordfilt2(CIM,sum(support(:)),support);

% determine the locations where the max over the neigh or 11x11
% corresponds to the cim values
[loc_x,loc_y]=find((cm==maxima).*(CIM>0));

% plot the corner features
figure;
imshow(im);
hold on;
plot(loc_y, loc_x, 'r*');
title('Harris Corner Detection');
