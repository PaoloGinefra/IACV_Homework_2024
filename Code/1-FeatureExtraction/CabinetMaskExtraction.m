close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Convert the image to HSV
im_hsv = rgb2hsv(im);

%% Plot all the channels
figure;
subplot(2, 3, 1);
imshow(im(:, :, 1));
title('Red Channel');

subplot(2, 3, 2);
imshow(im(:, :, 2));
title('Green Channel');

subplot(2, 3, 3);
imshow(im(:, :, 3));
title('Blue Channel');

subplot(2, 3, 4);
imshow(im_hsv(:, :, 1));
title('Hue Channel');

subplot(2, 3, 5);
imshow(im_hsv(:, :, 2));
title('Saturation Channel');

subplot(2, 3, 6);
imshow(im_hsv(:, :, 3));
title('Value Channel');
impixelinfo;

%% Plot the channels as 3d surfeces
figure;
subplot(2, 3, 1);
surf(im(:, :, 1), 'EdgeColor', 'none', 'FaceColor', 'interp');
title('Red Channel');

subplot(2, 3, 2);
surf(im(:, :, 2), 'EdgeColor', 'none', 'FaceColor', 'interp');
title('Green Channel');

subplot(2, 3, 3);
surf(im(:, :, 3), 'EdgeColor', 'none', 'FaceColor', 'interp');
title('Blue Channel');

subplot(2, 3, 4);
surf(im_hsv(:, :, 1), 'EdgeColor', 'none', 'FaceColor', 'interp');
title('Hue Channel');

subplot(2, 3, 5);
surf(im_hsv(:, :, 2), 'EdgeColor', 'none', 'FaceColor', 'interp');
title('Saturation Channel');

subplot(2, 3, 6);
surf(im_hsv(:, :, 3), 'EdgeColor', 'none', 'FaceColor', 'interp');
title('Value Channel');

%% Mask the saturation channel with a threshold
sat_threshold = 0.3;
sat_mask = im_hsv(:, :, 2) > sat_threshold;

green_threshold = 30;
green_mask = im(:, :, 2) > green_threshold;

sat_and_green_mask = sat_mask & green_mask;

erosion_kernel_size = 20;
erosion_kernel = strel('disk', erosion_kernel_size);
sat_and_green_mask_eroded = imerode(sat_and_green_mask, erosion_kernel);

dilation_kernel_size_increment = 20;
dilation_kernel = strel('disk', erosion_kernel_size + dilation_kernel_size_increment);
sat_and_green_mask_opened = imdilate(sat_and_green_mask_eroded, dilation_kernel);

gaussian_kernel_size = 50;
gaussian_sigma = 15;
sat_and_green_mask_opened_smoothed = imgaussfilt(double(sat_and_green_mask_opened), gaussian_sigma);

%% Plot the masks
figure;
nrows = 2;
ncols = 3;
subplot(nrows, ncols, 1);
imshow(sat_mask);
title(['Masked Saturation Channel - th:', num2str(sat_threshold)]);

subplot(nrows, ncols, 2);
imshow(green_mask);
title(['Masked Green Channel - th:', num2str(green_threshold)]);

subplot(nrows, ncols, 3);
imshow(sat_and_green_mask);
title('Saturation and Green Mask');

subplot(nrows, ncols, 4);
imshow(sat_and_green_mask_eroded);
title(['Eroded Mask - kernel size:', num2str(erosion_kernel_size)]);

subplot(nrows, ncols, 5);
imshow(sat_and_green_mask_opened);
title(['Opened Mask - kernel size:', num2str(erosion_kernel_size)]);

subplot(nrows, ncols, 6);
imshow(sat_and_green_mask_opened_smoothed);
title(['Opened Mask Smoothed - kernel size:', num2str(gaussian_kernel_size), ' sigma:', num2str(gaussian_sigma)]);

impixelinfo;

%% Mask the original image
im_masked = im;
im_masked(:, :, 1) = uint8(double(im_masked(:, :, 1)) .* sat_and_green_mask_opened_smoothed);
im_masked(:, :, 2) = uint8(double(im_masked(:, :, 2)) .* sat_and_green_mask_opened_smoothed);
im_masked(:, :, 3) = uint8(double(im_masked(:, :, 3)) .* sat_and_green_mask_opened_smoothed);

figure;
imshow(im_masked);
title('Masked Image');

%% Save the mask
image_mask = sat_and_green_mask_opened_smoothed;
save('./image_mask.mat', 'image_mask');
