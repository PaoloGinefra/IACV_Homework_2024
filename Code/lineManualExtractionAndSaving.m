% This file makes you draw all the necessary lines on the image and save them in a .mat file.
% The lines are saved in the following format:
% - [lineName]_points: a 2x2n matrix where n is the number of lines. Each pair of columns represent the start and end points of a line.
% - [lineName]_lines: a 3xn matrix where each column represents a line in homogeneous coordinates.

% The lines drawn are:
% - 4 x "h" lines
% - 4 x "m" lines
% - 3 x "l" lines

% To load the lines, you can use the following code:
% load('pathToMat/lines.mat');

close all;
clc;
clear;

%% Load the image
im = imread('../Assignment/Homework Image.jpg');

%% show the image
figure;
imshow(im);
title('Original Image');

%% get m lines using the draw line function
[h_points, h_lines] = drawLines('h', 4);
[m_points, m_lines] = drawLines('m', 4);
[l_points, l_lines] = drawLines('l', 3);

save('lines.mat', 'h_points', 'h_lines', 'm_points', 'm_lines', 'l_points', 'l_lines');
load(lines.mat);

function [points, lines] = drawLines(name, n_lines)
points = nan(2, 2 * n_lines);
lines = nan(3, n_lines);

for i = 1:n_lines
    figure(1);
    title(['Select the end points of the ' name ' lines']);
    segment = drawline("Color", "r");
    
    points(:, i:i+1) = segment.Position';
    a = [segment.Position(1,:)'; 1];
    b = [segment.Position(2,:)'; 1];
    l = cross(a, b);
    l = l / norm(l);
    lines(:, i) = l;
    
    midPoint = mean(segment.Position, 1);
    text(midPoint(1), midPoint(2), [name num2str(i)], 'Color', 'r', 'FontSize', 14);
end
end
