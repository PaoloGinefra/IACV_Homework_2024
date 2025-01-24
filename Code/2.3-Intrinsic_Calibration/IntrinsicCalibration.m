close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Load the lines for visualization purposes
load('../lines.mat');

%% load the line at infinity
load('../VanishingLineHorizontalPlane.mat');

%% Load the vanishing points
load('../vanishingPoints.mat');

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

l_v = l_intersection / l_intersection(3);
l_v_x = l_v(1);
l_v_y = l_v(2);
l_v_w = l_v(3);

m_v = m_intersection / m_intersection(3);
m_v_x = m_v(1);
m_v_y = m_v(2);
m_v_w = m_v(3);

v = h_intersection/h_intersection(3);
v_x = v(1);
v_y = v(2);
v_w = v(3);

H = inv(H_metric);
h1 = H(:, 1) / H(3, 1);
h2 = H(:, 2) / H(3, 2);
h1_x = h1(1);
h1_y = h1(2);
h1_w = h1(3);
h2_x = h2(1);
h2_y = h2(2);
h2_w = h2(3);

%% Plot the image
figure;
imshow(im);
hold on;
for i = 1:size(h_lines, 2)
    line(h_points(1, 2*i-1:2*i), h_points(2, 2*i-1:2*i), 'Color', 'b', 'LineWidth', 2);
end


% plot h1 and h2
plot([h1_x, h2_x], [h1_y, h2_y], 'r', 'LineWidth', 2);
% plot the vanishing point
plot(v_x, v_y, 'bx', 'MarkerSize', 10, 'LineWidth', 2);
%plot m_v and l_v
plot(l_v_x, l_v_y, 'gx', 'MarkerSize', 10, 'LineWidth', 2);
plot(m_v_x, m_v_y, 'gx', 'MarkerSize', 10, 'LineWidth', 2);

l_h_inf_prime = l_h_inf_prime / l_h_inf_prime(3);
l_h_inf_prime_x = l_h_inf_prime(1);
l_h_inf_prime_y = l_h_inf_prime(2);
l_h_inf_prime_w = l_h_inf_prime(3);

% A = [
%     h1_x * h2_x, h1_x, h1_y, 1;
%     h1_x^2 - h2_x^2, h1_x-h2_x, h1_y - h2_y, 0;
%     v_x, 1, 0, 0;
%     0, 0, 1, 0
%     0, 0, 0, 1
%     ];

% b = [
%     -h1_y*h2_y;
%     -h1_y^2 + h2_y^2;
%     l_h_inf_prime_x;
%     l_h_inf_prime_y - v_y;
%     1
%     ];

% A = [
%     h1_x * h2_x, h1_x, h1_y*h2_y, h1_y;
%     h1_x^2 - h2_x^2, h1_x-h2_x, h1_y^2 - h2_y^2, h1_y - h2_y;
%     v_x, 1, 0, 0;
%     0, 0, v_y, 1
%     ];

% b = [
%     -1;
%     0;
%     l_h_inf_prime_x;
%     l_h_inf_prime_y;
%     ];

% A = [
%     h1_x*h2_x, h1_x*h2_w, h1_w*h2_w;
%     h1_x^2 - h2_x^2, h1_x*h1_w - h2_x*h2_w, h1_w^2 - h2_w^2;
%     v_x, v_w, 0
%     0, 0, v_w
%     ];

% b = [
%     -1;
%     0;
%     l_h_inf_prime_x;
%     l_h_inf_prime_w;
%     ];

% A = [
%     h1_x * h2_x, h1_x, h1_y*h2_y, h1_y;
%     h1_x^2 - h2_x^2, h1_x-h2_x, h1_y^2 - h2_y^2, h1_y - h2_y;
%     v_x * h1_x, v_x, v_y * h1_y, v_y;
%     v_x * h2_x, v_x, v_y * h2_y, v_y
%     ];

% b = [
%     -1;
%     0;
%     -1;
%     -1;
%     ];

A = [
    h1_x * h2_x, h1_x + h2_x, h1_y+h2_y, 1;
    h1_x^2 - h2_x^2, 2 *(h1_x-h2_x), 2*(h1_y - h2_y), 0;
    v_x * h1_x, v_x+h1_x, v_y+h1_y, 1;
    v_x * h2_x, v_x+h2_x, v_y+h2_y, 1;
    ];

b = [
    -h1_y*h2_y;
    -h1_y^2 + h2_y^2;
    -v_y*h1_y;
    -v_y*h2_y;
    ];

x = pinv(A) * b;

% omega = [
%     x(1), 0, x(2);
%     0, 1, x(3);
%     0, 0, x(4)
%     ];

% omega = [
%     x(1), 0, x(2);
%     0, x(3), x(4);
%     0, 0, 1
%     ];

% omega = [
%     x(1), 0, x(2);
%     0, 1, 0;
%     0, 0, x(3)
%     ];

% omega = [
%     x(1), 0, x(2);
%     0, x(3), x(4);
%     0, 0, 1
%     ];

omega = [
    x(1), 0, x(2);
    0, 1, x(3);
    x(2), x(3), x(4)
    ];
% omega= omega/omega(3,3);

% omega = omega/x(3);

disp('h1T Omega h2');
disp(h1'*omega*h2);

disp('h1T Omega h1 - h2T Omega h2');
disp(h1'*omega*h1 - h2'*omega*h2);

disp('vT Omega h1');
disp(v'*omega*h1);

disp('vT Omega h2');
disp(v'*omega*h2);

disp('Eigs of Omega');
disp(eig(omega));

%% Compute K using the cholesky decomposition of omega inverse

alfa = sqrt(omega(1,1));
u0 = -omega(1,3)/(alfa^2);
v0 = -omega(2,3);
fy = sqrt(omega(3,3) - (alfa^2)*(u0^2) - (v0^2));
fx = fy / alfa;

% build K using the parametrization
K = [fx 0 u0; 0 fy v0; 0 0 1];

disp('Sanity Check [K]: omega^-1 = KK^T');
omega_inv = inv(omega);
omega_inv = omega_inv/omega_inv(3,3);
omega_inv_K = K*K';
disp(norm(omega_inv - omega_inv_K));
% omega_inv = inv(omega);

% K = chol(omega_inv, 'lower');

disp('K');
disp(K);

%% Save the intrinsic matrix
save('../K_metric.mat', 'K');

rs = K\H;
disp('Rs');
disp(rs)

disp('Sanity Check [K]: Orthogonality of world coordinate axes');
disp(rs(:,1)'*rs(:,2)/(norm(rs(:,1))*norm(rs(:,2))));

%% make a 3d plot of colums of rs
figure;
plot3(rs(1, 3), rs(2, 3), rs(3, 3),  '+', 'LineWidth', 2);
grid on;

% draw two arrows starting from the point and in the direction of the first 2 cols of rs
hold on;
quiver3(rs(1, 3), rs(2, 3), rs(3, 3), rs(1, 1), rs(2, 1), rs(3, 1), 0, 'LineWidth', 2);
quiver3(rs(1, 3), rs(2, 3), rs(3, 3), rs(1, 2), rs(2, 2), rs(3, 2), 0, 'LineWidth', 2);

%Plot the plane of the span of the first 2 columns of rs
% create a grid
[X, Y] = meshgrid(0:1:l1_length, 0:1:0.35*l1_length);

grid_points = [X(:), Y(:)];
grid_points_cart = rs(:, 1:2) * grid_points' + rs(:, 3);
X = reshape(grid_points_cart(1, :), size(X));
Y = reshape(grid_points_cart(2, :), size(Y));
Z = reshape(grid_points_cart(3, :), size(X));
surf(X, Y, Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

%% Plot the lines
world_l_points = rs(:, 1:2) * l_points_metric + rs(:, 3);
plot3(world_l_points(1, 1:2), world_l_points(2, 1:2), world_l_points(3, 1:2), 'r', 'LineWidth', 2);

world_m_points = rs(:, 1:2) * m_points_metric + rs(:, 3);
plot3(world_m_points(1, 1:2), world_m_points(2, 1:2), world_m_points(3, 1:2), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 3:4), world_m_points(2, 3:4), world_m_points(3, 3:4), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 5:6), world_m_points(2, 5:6), world_m_points(3, 5:6), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 7:8), world_m_points(2, 7:8), world_m_points(3, 7:8), 'g', 'LineWidth', 2);
% axis of same scale
axis equal;

% add axis labels
xlabel('X');
ylabel('Y');
zlabel('Z');


%% Compute the height h
height = (1/l1_length - 1/l2_length) * vecnorm(rs(:, 3));
disp('Height');
disp(height);