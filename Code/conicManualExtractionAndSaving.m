close all;
clc;
clear;

%% Load the image
im = imread('../Assignment/Homework Image.jpg');

%% show the image
figure;
imshow(im);
title('Select at least 5 points of the desired conic');
hold on;

%% get points on the conic using getpts
[xs, ys] = getpts(figure(1));
n_points = length(xs);

if(n_points < 5)
    error('You need to select at least 5 points to define a conic');
end

%% Construct the A matrix
A = zeros(n_points, 6);
A(:, 1) = xs.^2;
A(:, 2) = xs .* ys;
A(:, 3) = ys.^2;
A(:, 4) = xs;
A(:, 5) = ys;
A(:, 6) = ones(n_points, 1);

%% Solve the system
[~, ~, V] = svd(A);
conic = V(:, end);

A * conic

%% Construct the conic matrix
C = [conic(1), conic(2) / 2, conic(4) / 2;
    conic(2) / 2, conic(3), conic(5) / 2;
    conic(4) / 2, conic(5) / 2, conic(6)];

%% Plot the conic
conicPlot = zeros(size(im, 1), size(im, 2));

for i = 1:size(im, 1)
    for j = 1:size(im, 2)
        point = [j, i, 1];
        conicPlot(i, j) = point * C * point';
    end
end

% plot the conic
contour(conicPlot, [0, 0], 'r', 'LineWidth', 2);

%plot the points
plot(xs, ys, 'rx', 'MarkerSize', 10, 'LineWidth', 2);

title('Conic');

%% Save the conic
conicName = input('Enter the name of the conic: ', 's');
save([conicName, '.mat'], 'C');

