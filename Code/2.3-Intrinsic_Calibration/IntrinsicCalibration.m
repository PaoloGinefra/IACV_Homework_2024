close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Load the lines for visualization purposes
load('../lines.mat');

%% load the line at infinity
load('../VanishingLineHorizontalPlane.mat');

%% load conic
% load('./CircleC.mat');
load('../1-FeatureExtraction/ConicExtraction/ExtractedConic.mat');

%% Load the metric rectification homography
load('../H_metric.mat');

%% Compute the vanishing point ortogonal to the rectified face
%%Find the intersection of the h lines
[~, ~, V] = svd(h_lines');
h_intersection = V(:, end);
h_intersection = h_intersection / h_intersection(3);
h_intersection_euclidian = h_intersection(1:2);

%% plot the h lines
figure;
imshow(im);
hold on;
for i = 1:size(h_lines, 2)
    line(h_points(1, 2*i-1:2*i), h_points(2, 2*i-1:2*i), 'Color', 'b', 'LineWidth', 2);
end

%% plot the intersection point
plot(h_intersection_euclidian(1), h_intersection_euclidian(2), 'bx', 'MarkerSize', 10, 'LineWidth', 2);



v = h_intersection;
v_x = v(1);
v_y = v(2);

H = inv(H_metric);
h1 = H(:, 1) / H(3, 1);
h2 = H(:, 2) / H(3, 2);
h1_x = h1(1);
h1_y = h1(2);
h2_x = h2(1);
h2_y = h2(2);

l_h_inf_prime = l_h_inf_prime / l_h_inf_prime(3);
l_h_inf_prime_x = l_h_inf_prime(1);
l_h_inf_prime_y = l_h_inf_prime(2);

A = [
    h1_x * h2_x, h1_x, h1_y, 1;
    h1_x^2 - h2_x^2, h1_x-h2_x, h1_y - h2_y, 0;
    v_x, 1, 0, 0;
    0, 0, 1, 0
    0, 0, 0, 1
    ];

b = [
    -h1_y*h2_y;
    -h1_y^2 + h2_y^2;
    l_h_inf_prime_x;
    l_h_inf_prime_y - v_y;
    1
    ];

x = A\b;

omega = [
    x(1), 0, x(2);
    0, 1, x(3);
    0, 0, x(4)
    ];

%omega = omega / x(3);

%% Compute K using the cholesky decomposition of omega inverse
omega_inv = inv(omega);
omega_inv = omega_inv / omega_inv(3, 3);

K = chol(omega_inv, 'lower');

disp('K');
disp(K);

%% Save the intrinsic matrix
save('../K_metric.mat', 'K');