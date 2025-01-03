close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Load the edge mask
load('./edges_2025-01-03_15-08-15.mat');
load('./edges_inner_2025-01-03_16-14-00.mat');

im_edges = im_edges + inner_edges;

%% Load the corner points
load('./corner_points_2025-01-03_15-07-18.mat');

%% Extract the coordinates of the white pixel in the edge image
[x, y] = find(im_edges > 0.5);
edge_points = [y, x];

%% Plot the edge points
figure;
imshow(im);
alpha(0.5);
hold on;
plot(edge_points(:, 1), edge_points(:, 2), 'r.', 'MarkerSize', 1);
plot(corner_points(:, 1), corner_points(:, 2), 'g.', 'MarkerSize', 10);
hold off;
title('Edge Points');

%% RANSAC parameters
numIterations = 50000;
inlierThreshold = 5;
numPoints = 2;
eps = 3;
minLength = 150;
lambda_length = 0;%0.0003;


nModels = 20;
bestModels = zeros(nModels, 2, 2);

for m = 1:nModels
    bestModel = [];
    bestScore = -Inf;
    for i = 1:numIterations
        %% Randomly select two points
        idx = randperm(size(corner_points, 1), numPoints);
        p1 = corner_points(idx(1), :);
        p2 = corner_points(idx(2), :);
        
        %% Compute the score
        [score, nInliners] = computeScore(p1, p2, edge_points, eps, minLength, lambda_length);
        
        if score > bestScore
            bestScore = score;
            bestModel = [p1; p2];
            disp(['Iteration: ', num2str(i), ' Score: ', num2str(score), ' nInliners: ', num2str(nInliners), ' P1 - ', num2str(p1), ' P2 - ', num2str(p2)]);
            % Plot the best model
            figure(34);
            imshow(im);
            hold on;
            % Plot all the best models
            for j = 1:m
                p1_ = bestModels(j, 1, :);
                p2_ = bestModels(j, 2, :);
                p1_ = reshape(p1_, 1, 2);
                p2_ = reshape(p2_, 1, 2);
                line([p1_(1), p2_(1)], [p1_(2), p2_(2)], 'Color', 'b', 'LineWidth', 2);
            end
            drawModel(p1, p2, eps);
            plot(edge_points(:, 1), edge_points(:, 2), 'r.', 'MarkerSize', 2);
            plot(corner_points(:, 1), corner_points(:, 2), 'g.', 'MarkerSize', 10);
            hold off;
        end
    end
    bestModels(m, :, :) = bestModel;
    
    % Remove the model from the corner points
    corner_points = setdiff(corner_points, bestModel, 'rows');
end


function [score, nInliners] = computeScore(p1, p2, points, eps, minLength, lambda_length)
dir = p2 - p1;
l = norm(dir);

if(l < minLength)
    score = -Inf;
    nInliners = 0;
    return;
end

dir = dir / l;

projectionMatrix = [dir; -dir(2), dir(1)];
deltas = points - p1;
projections = projectionMatrix * deltas';


projections = projections(:, projections(1, :) > 0 & projections(1, :) < l);
projections = projections(:, projections(2, :) > -eps & projections(2, :) < eps);

nInliners = size(projections, 2);

%% Compute the score
projections(2, :) = -abs(projections(2, :)) + eps;
projections(2, :) = max(projections(2, :), 0);
projections(2, :) = projections(2, :)/eps;

score = sum(projections(2, :)) + lambda_length * l;
end

function drawModel(p1, p2, eps)
dir = p2 - p1;
l = norm(dir);
dir = dir / l;

line([p1(1), p2(1)], [p1(2), p2(2)], 'Color', 'r', 'LineWidth', eps);
alpha(0.5);
plot([p1(1) - eps * dir(2), p2(1) - eps * dir(2)], [p1(2) + eps * dir(1), p2(2) + eps * dir(1)], 'b--', 'LineWidth', 1);
plot([p1(1) + eps * dir(2), p2(1) + eps * dir(2)], [p1(2) - eps * dir(1), p2(2) - eps * dir(1)], 'b--', 'LineWidth', 1);
end