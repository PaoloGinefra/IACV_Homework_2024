close all;
clc;
clear;

%% Load Image
im = imread('../../Assignment/Homework Image.jpg');

%% Convert to grayscale
im = im(:, :, 2);
im = double(im)/255;

figure;
imshow(im);
title('Original Image');

%% SVD decomposition
[U, S, V] = svd(im);

%% Plot the singular values
figure;
subplot(1, 3, 1);
semilogy(diag(S));

% Plot the cumulative sum of singular values
subplot(1, 3, 2);
plot(cumsum(diag(S))/sum(diag(S)));

% Plot the explained variance
subplot(1, 3, 3);
plot(cumsum(power(diag(S), 2))/sum(power(diag(S), 2)));


%% Low rank approximation
figure;
ranks = [1, 5, 10, 20, 50, 100, 200, 500];
nRows = 2;
nCols = 4;

for i = 1:length(ranks)
    rank = ranks(i);
    im_approx = ComputeLowRankApproximation(U, S, V, rank);
    
    subplot(nRows, nCols, i);
    imshow(im_approx);
    title(['Rank: ', num2str(rank)]);
end

%% Rank-1 component
figure;
for i = 1:length(ranks)
    rank1Component = ComputeRank1Component(U, S, V, i);
    
    subplot(nRows, nCols, i);
    imshow(rank1Component);
    title(['Rank-1 Component: ', num2str(i)]);
end


function [im_approx] = ComputeLowRankApproximation(U, S, V, rank)
% Low rank approximation
im_approx = U(:, 1:rank) * S(1:rank, 1:rank) * V(:, 1:rank)';
end

function [rank1Component] = ComputeRank1Component(U, S, V, n)
% Compute the rank-1 component
rank1Component = U(:, n) * S(n, n) * V(:, n)';
end
