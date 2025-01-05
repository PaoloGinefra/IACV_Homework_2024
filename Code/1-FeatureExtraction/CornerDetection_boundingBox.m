close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

figure(47);
imshow(im);
impixelinfo;
title('Original Image');

%% Convert the image to HSV
im_hsv = rgb2hsv(im);

%% Convert the image to grayscale
im_gray = im_hsv(:, :, 2);%%rgb2gray(im);

%% Rescale the image
% min_val = 79;
% max_val = 120;
% im_gray = imadjust(im_gray, [min_val/255, max_val/255], [0, 1]);

%% Apply Cannys edge detection
% check type of image
if ~isa(im_gray, 'double')
    y=double(im_gray) / 255.0;
else
    y = im_gray;
end
figure(1), imshow(y),title('original image');

dx = [-1 0 1; -1 0 1; -1 0 1];   % Derivative masks
dy = dx';

Ix = conv2(y, dx, 'same');      % Image derivatives
Iy = conv2(y, dy, 'same');

Gradients = sqrt(Ix.^2 + Iy.^2);
figure(2), imshow(Gradients, []), title('Gradients');

% set the parameter for Gaussian convolution used in Harris Corner Detector
SIGMA_gaussian=5
g = fspecial('gaussian',max(1,fix(3*SIGMA_gaussian)+1), SIGMA_gaussian);

Ix2 = conv2(Ix.^2, g, 'same'); % Smoothed squared image derivatives
Iy2 = conv2(Iy.^2, g, 'same');
Ixy = conv2(Ix.*Iy, g, 'same');

% cim = det(M) - k trace(M)
k = 0.04;
detM = Ix2.*Iy2 - Ixy.^2;
traceM = Ix2 + Iy2;
cim = detM - k * traceM;
%cim = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps);

figure(); imshow(detM, []), title('det(M)');
figure(); imshow(traceM, []), title('trace(M)');

%% Crop the border of the image
BORDER=20;
cim(1:BORDER,:)=0;
cim(end-BORDER:end,:)=0;
cim(:,end-BORDER:end)=0;
cim(:,1:BORDER)=0;


%% Thresholding the cim
T=mean(cim(:));
CIM=cim;
CIM(find(cim<T))=0;
% similarly one could use the Otzu method

figure(), imshow(CIM,[]),title('Harris measure');
colorbar

figure(34), mesh(CIM),title('Harris measure');
colorbar

%% perform nonlocal maximum suppression on the thresholded measure
% this value needs to be adjsted also depending on the image size
support=true(11);
% compute maximum over a square neighbor of size 11 x 11
maxima=ordfilt2(CIM,sum(support(:)),support);
% determine the locations where the max over the neigh or 11 x 11 corresponds to the cim values
[loc_x,loc_y]=find((cim==maxima).*(CIM>0));
indx = find((cim==maxima).*(CIM>0));

disp('local maxima on the cim measure displayed as a mesh')
figure(34),
hold on
plot3(loc_y,loc_x, cim(indx), 'g+', 'LineWidth', 2)
hold off
view(gca,[-66.3 42.8]);


% draw a cross on the image in the local maxima
figure(), imshow(y,[]), hold on,
plot(loc_y,loc_x,'g+', 'LineWidth', 4)
title('Local maxima of Harris measure')

%% Save the corner points
corner_points = [loc_y, loc_x];
currentTime = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
save(['./corner_points_', currentTime, '.mat'], 'corner_points');