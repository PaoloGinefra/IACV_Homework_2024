%% This file calculates the vanishing line of the horizontal plane using the lines extracted from the image.
% The lines are loaded from a file in the following format:
% - [lineName]_points: a 2x2n matrix where n is the number of lines. Each pair of columns represent the start and end points of a line.
% - [lineName]_lines: a 3xn matrix where each column represents a line in homogeneous coordinates.
%
% The lines used are:
% - 3 x "l" lines
% - 4 x "m" lines
%
% The vanishing line is saved in a .mat file in the following format:
% - l_h_inf_prime: the vanishing line of the horizontal plane in homogeneous coordinates.
%
%This script also saves the vanishing points of the lines in a .mat file in the following format:
% - l_intersection: the vanishing point of the "l" lines in homogeneous coordinates.
% - m_intersection: the vanishing point of the "m" lines in homogeneous coordinates.
%
% To load the vanishing line, you can use the following code:
% load('pathToMat/vanishingLineHorizontalPlane.mat');
%
% To load the vanishing points, you can use the following code:
% load('pathToMat/vanishingPoints.mat');


close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Load the lines
load('../2.0-ManualLineExtraction/lines.mat');

%%Find the intersection of the l lines
[~, ~, V] = svd(l_lines');
l_intersection = V(:, end);
l_intersection_euclidian = l_intersection(1:2) / l_intersection(3);

%%Find the intersection of the m lines
[~, ~, V] = svd(m_lines');
m_intersection = V(:, end);
m_intersection_euclidian = m_intersection(1:2) / m_intersection(3);

%%Plot the j and m lines
figure;
imshow(im);
hold on;
for i = 1:size(l_lines, 2)
    line(l_points(1, 2*i-1:2*i), l_points(2, 2*i-1:2*i), 'Color', 'r', 'LineWidth', 2);
end

for i = 1:size(m_lines, 2)
    line(m_points(1, 2*i-1:2*i), m_points(2, 2*i-1:2*i), 'Color', 'g', 'LineWidth', 2);
end

%%Plot the intersection points
plot(l_intersection_euclidian(1), l_intersection_euclidian(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(m_intersection_euclidian(1), m_intersection_euclidian(2), 'gx', 'MarkerSize', 10, 'LineWidth', 2);

%%Calculate the vanishing line
l_h_inf_prime = cross(l_intersection, m_intersection);

%%Plot the vanishing line
x_min =  min(1, min(m_intersection_euclidian(1), l_intersection_euclidian(1)));
x_max = max(size(im, 2), max(m_intersection_euclidian(2), l_intersection_euclidian(2)));
x = x_min:x_max;
y = (-l_h_inf_prime(3) - l_h_inf_prime(1) * x) / l_h_inf_prime(2);
plot(x, y, 'b--', 'LineWidth', 2);
text(mean(x), mean(y), 'Vanishing Line', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold');
title('Vanishing Line');

%% Save the vanishing line
save('vanishingLineHorizontalPlane.mat', 'l_h_inf_prime');

%% Save the vanishing points
save('vanishingPoints.mat', 'l_intersection', 'm_intersection');