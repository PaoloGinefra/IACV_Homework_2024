close all;
clc;
clear;

%% Load the image
im = imread('../../../Assignment/Homework Image.jpg');

%% Cut out the region of interest
% Define the region of interest
BL_corner = [328, 570];
UR_corner = [773, 505];

% Cut out the region of interest
im = im(UR_corner(2):BL_corner(2), BL_corner(1):UR_corner(1), :);
figure;
imshow(im);
title('Cropped Image');

% Save the plot for the report
saveas(gcf, 'CroppedImage.png');

% convert to hsv
im_hsv = rgb2hsv(im);

% convert to grayscale
im_gray = im_hsv(:, :, 3);

figure;
imshow(im_gray);
title('Grayscale Image');

% Save the plot for the report
saveas(gcf, 'GrayscaleImage.png');

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

% apply Cannys edge detection
threshold = [0, 0.15];
sigma = 2;
im_edges = edge(im_gray, 'canny', threshold, sigma);

% Display the results
figure;
subplot(2, 1, 1);
imshow(im_gray);
title('Grayscale Image');

subplot(2, 1, 2);
imshow(im_edges);
title(['Canny Edges - Threshold: ', num2str(threshold), ' Sigma: ', num2str(sigma)]);

% Save the plot for the report
saveas(gcf, 'CannyEdges.png');

% Blur the edges
sigma = 10;
im_edges_blurred = imgaussfilt(double(im_edges), sigma);

% Display the results
figure;
subplot(3, 1, 1);
imagesc(im_edges_blurred);
title(['Blurred Edges - Sigma: ', num2str(sigma)]);
impixelinfo;

% threshold the blurred edges
threshold = 0.1;
im_edges_thresholded = im_edges_blurred > threshold;

% Display the results
subplot(3, 1, 2);
imagesc(im_edges_thresholded);
title(['Thresholded Edges - Threshold: ', num2str(threshold)]);
impixelinfo;

%mask the edges
im_edges = im_edges .* (1-im_edges_thresholded);

% Display the results
subplot(3, 1, 3);
imshow(im_edges);
title('Masked Edges');

% Save the plot for the report
saveas(gcf, 'MaskedEdges.png');


%% RANSAC to find the conic
% points from edge detection
[y, x] = find(im_edges(1:60, :));
points = [x'; y'];

% RANSAC parameters
n = 10000;
t = 0.01; % inlier threshold
d = 0; % percentage of inliers
s = 5; % number of points to fit the conic

% Visualization setup
figure;
subplot(3, 1, 1);
imshow(im);

subplot(3, 1, 2);
imshow(im_edges);

subplot(3, 1, 3);
hold on;

bestConic = zeros(3, 3);
bestScore = 0;

% RANSAC
for step = 1:n
    % selct 5 random points
    idx = randperm(size(points, 2), s);
    sample = points(:, idx);
    
    % fit a conic
    C = conicFit(sample);
    
    % calculate the distance
    dist = conicDistance(C, points);
    
    % find the inliers
    inliers = find(dist < t);
    
    % check if the number of inliers is greater than the threshold
    if length(inliers) > d * length(points)
        % fit the conic to all the inliers
        inlierPoints = points(:, inliers);
        C = conicFit(inlierPoints);
        
        % calculate the score
        score = length(inliers);
        
        % update the best conic
        if score > bestScore
            bestScore = score;
            bestConic = C;
            
            % plot the conic
            conicPlot = zeros(size(im, 1), size(im, 2));
            
            for i = 1:size(im, 1)
                for j = 1:size(im, 2)
                    point = [j, i, 1];
                    conicPlot(i, j) = point * C * point';
                end
            end
            
            subplot(3, 1, 3);
            % plot the image
            imshow(im);
            
            % plot the conic
            contour(conicPlot, [0, 0], 'r', 'LineWidth', 2);
            
            % plot the points
            plot(inlierPoints(1, :), inlierPoints(2, :), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
            
            %plot the orignal sample points
            plot(sample(1, :), sample(2, :), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
            
            title(['Conic at iteration: ', num2str(step), ' with score: ', num2str(score)]);
            
            pause(0.1);
        end
    end
    
    
end

% save the plot
saveas(gcf, 'ConicRANSAC.png');


% plot the final conic
conicPlot = zeros(size(im, 1), size(im, 2));

for i = 1:size(im, 1)
    for j = 1:size(im, 2)
        point = [j, i, 1];
        conicPlot(i, j) = point * bestConic * point';
    end
end

figure;
imshow(im);
hold on;
contour(conicPlot, [0, 0], 'r', 'LineWidth', 2);
title('Final Conic');

% Save the plot for the report
saveas(gcf, 'FinalConic.png');

% save the conic
conicName = input('Enter the name of the conic: ', 's');
%Adjust the conic for the initial cropping
translation = [1, 0, -BL_corner(1); 0, 1, -UR_corner(2); 0, 0, 1];
C = translation' * bestConic * translation;
save([conicName, '.mat'], 'C');

function [C] =  conicFit(points)
if size(points, 2) < 5
    error('You need to select at least 5 points to define a conic');
end

% Construct the A matrix
A = zeros(size(points, 2), 6);
xs = points(1, :);
ys = points(2, :);

A(:, 1) = xs.^2;
A(:, 2) = xs .* ys;
A(:, 3) = ys.^2;
A(:, 4) = xs;
A(:, 5) = ys;
A(:, 6) = ones(size(points, 2), 1);

% Solve the system
[~, ~, V] = svd(A);
conic = V(:, end);

% Construct the conic matrix
C = [conic(1), conic(2) / 2, conic(4) / 2;
    conic(2) / 2, conic(3), conic(5) / 2;
    conic(4) / 2, conic(5) / 2, conic(6)];
end

function [dist] = conicDistance(C, points)
homogeneousPoints = [points; ones(1, size(points, 2))];
dist = zeros(1, size(points, 2));

for i = 1:size(points, 2)
    dist(i) = homogeneousPoints(:, i)' * C * homogeneousPoints(:, i);
end

dist = abs(dist);
end
