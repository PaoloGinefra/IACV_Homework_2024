close all;
clc;
clear;

%% Load the image
im = imread('../../Assignment/Homework Image.jpg');

%% Load the lines for visualization purposes
load('../2.0-ManualLineExtraction/lines.mat');

%% load the line at infinity
load('../2.1-Vanishing_Line_Extraction/VanishingLineHorizontalPlane.mat');

%% Load the vanishing points
load('../2.1-Vanishing_Line_Extraction/vanishingPoints.mat');

%% load conic
% load('./CircleC.mat');
load('../1-FeatureExtraction/ConicExtraction/ExtractedConic.mat');

%% Load the metric rectification homography
load('../2.2-Metric_Rectification/H_metric.mat');

%% Compute the vanishing point ortogonal to the rectified face
%% Find the intersection of the h lines
[~, ~, V] = svd(h_lines');
h_intersection = V(:, end);
h_intersection = h_intersection / h_intersection(3);

v = h_intersection/h_intersection(3);

H = inv(H_metric);

%% Compute the omega matrix
omega = get_IAC(l_h_inf_prime, v, h_intersection, m_intersection, H);

%% Sanity checks
assert(norm(omega - omega') < 1e-10);
disp('Sanity check[Omega] - Omega is symmetric: ✅');

h1 = H(:, 1);
h2 = H(:, 2);

assert(abs(h1'*omega*h2) < 1e-10);
disp('Sanity check[Omega] - h1T Omega h2 = 0: ✅');

assert(abs(h1'*omega*h1 - h2'*omega*h2) < 1e-10);
disp('Sanity check[Omega] - h1T Omega h1 = h2T Omega h2: ✅');

% disp(v'*omega*h1);
assert(abs(v'*omega*h1) < 1e-10);
disp('Sanity check[Omega] - vT Omega h1 = 0: ✅');

% disp(v'*omega*h2);
assert(abs(v'*omega*h2) < 1e-10);
disp('Sanity check[Omega] - vT Omega h2 = 0: ✅');

assert(all(eigs(omega) > 0));
disp('Sanity check[Omega] - Omega is positive definite: ✅');

%% Compute K from omega
alfa = sqrt(omega(1,1));
u0 = -omega(1,3)/(alfa^2);
v0 = -omega(2,3);
fy = sqrt(omega(3,3) - (alfa^2)*(u0^2) - (v0^2));
fx = fy / alfa;

% build K using the parametrization
K = [fx 0 u0; 0 fy v0; 0 0 1];

%% Sanity checks
omega_inv = inv(omega);
omega_inv = omega_inv/omega_inv(3,3);
omega_inv_K = K*K';
assert(norm(omega_inv - omega_inv_K) < 1e-9);
disp('Sanity Check [K]: omega^-1 = KK^T: ✅');

%% Sanity check on K: compute the basis and offset of the rectified plane and check the orthogonality of the basis
% Compute the basis and offset of the rectified plane
rs = K\H;
% Scale the basis so as to use the metric coordinates of the plane containg l1 (The upper face of the parallelepiped)
rs = rs / norm(rs(:, 1)) / l1_length;
r3 = cross(rs(:, 1), rs(:, 2));
r3 = r3 / norm(r3);

assert(abs(rs(:,1)'*rs(:,2)/(norm(rs(:,1))*norm(rs(:,2)))) < 1e-10);
disp('Sanity Check [K]: Orthogonality of the basis of the rectified plane (r1 ⊥ r2): ✅');

%% Print K
disp('K:');
disp(K);

%% Save the intrinsic matrix
save('./K.mat', 'K', 'rs', 'r3');

%% make a 3d plot of colums of rs
figure;
plot3(rs(1, 3), rs(2, 3), rs(3, 3),  '+', 'LineWidth', 2);
grid on;

% draw two arrows starting from the point and in the direction of the first 2 cols of rs
hold on;
quiver3(rs(1, 3), rs(2, 3), rs(3, 3), rs(1, 1), rs(2, 1), rs(3, 1), 0, 'LineWidth', 2);
quiver3(rs(1, 3), rs(2, 3), rs(3, 3), rs(1, 2), rs(2, 2), rs(3, 2), 0, 'LineWidth', 2);

%Plot the plane of the span of the first 2 columns of rs
%create a grid
[X, Y] = meshgrid(0:l1_length / 100:l1_length, 0:l1_length/100:0.35 * l1_length);

grid_points = [X(:), Y(:)];
grid_points_cart = rs(:, 1:2)* grid_points' + rs(:, 3);
X = reshape(grid_points_cart(1, :), size(X));
Y = reshape(grid_points_cart(2, :), size(Y));
Z = reshape(grid_points_cart(3, :), size(X));
surf(X, Y, Z, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

%% Plot the lines
world_l_points = rs(:, 1:2) * l_points_metric + rs(:, 3);
world_l_points(:, 3:6) = world_l_points(:, 3:6) * l1_length/l2_length;
plot3(world_l_points(1, 1:2), world_l_points(2, 1:2), world_l_points(3, 1:2), 'r', 'LineWidth', 2);
plot3(world_l_points(1, 3:4), world_l_points(2, 3:4), world_l_points(3, 3:4), 'r', 'LineWidth', 2);
plot3(world_l_points(1, 5:6), world_l_points(2, 5:6), world_l_points(3, 5:6), 'r', 'LineWidth', 2);

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

% plot a camera at the origin
plotCamera('Location', [0, 0, -0.2], 'Orientation', eye(3), 'Size', 0.1, 'Color', 'b');


%% Compute the height h
height = (1 - l1_length/l2_length) * rs(:, 3)' * r3;
disp('Height');
disp(height);

%% Plot the S curve
load('../S_points.mat');
S_points_world = rs(:, 1:2) * S_points' + rs(:, 3);
scaling = 1 - height / (2 * rs(:, 3)' * r3);

S_points_world = S_points_world * scaling;

plot3(S_points_world(1, :), S_points_world(2, :), S_points_world(3, :), 'b', 'LineWidth', 2);

%% Move everything to the origin
origin = world_l_points(:, 3);
i_prime = rs(:, 1) / norm(rs(:, 1));
j_prime = rs(:, 2) / norm(rs(:, 2));
k_prime = r3;

R_inv = [i_prime, j_prime, k_prime, origin;
    0, 0, 0, 1];

R = inv(R_inv);

world_l_points = R * [world_l_points; ones(1, size(world_l_points, 2))];
world_m_points = R * [world_m_points; ones(1, size(world_m_points, 2))];
S_points_world = R * [S_points_world; ones(1, size(S_points_world, 2))];
camera_pos = R * [0; 0; 0; 1];
camera_direction = R * [0; 0; 1; 1];
camera_direction = camera_direction / camera_direction(4);
camera_direction = camera_direction - camera_pos;
camera_direction = camera_direction / norm(camera_direction);
camera_direction = camera_direction(1:3);
camera_direction = camera_direction / norm(camera_direction);
%Compute camera orientation from camera direction for the plot camera
i_camera_direction = [camera_direction(2); -camera_direction(1); 0];
i_camera_direction = i_camera_direction / norm(i_camera_direction);

j_camera_direction = cross(camera_direction, i_camera_direction);
cameraOrientation = inv([i_camera_direction, j_camera_direction, camera_direction]);


figure;
plot3(world_l_points(1, 1:2), world_l_points(2, 1:2), world_l_points(3, 1:2), 'r', 'LineWidth', 2);
hold on;
plot3(world_l_points(1, 3:4), world_l_points(2, 3:4), world_l_points(3, 3:4), 'r', 'LineWidth', 2);
plot3(world_l_points(1, 5:6), world_l_points(2, 5:6), world_l_points(3, 5:6), 'r', 'LineWidth', 2);

plot3(world_m_points(1, 1:2), world_m_points(2, 1:2), world_m_points(3, 1:2), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 3:4), world_m_points(2, 3:4), world_m_points(3, 3:4), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 5:6), world_m_points(2, 5:6), world_m_points(3, 5:6), 'g', 'LineWidth', 2);
plot3(world_m_points(1, 7:8), world_m_points(2, 7:8), world_m_points(3, 7:8), 'g', 'LineWidth', 2);

plot3(S_points_world(1, :), S_points_world(2, :), S_points_world(3, :), 'b', 'LineWidth', 2);

% axis of same scale
axis equal;

% grid on
grid on;

% add axis labels
xlabel('X');
ylabel('Y');
zlabel('Z');

% plot camera noraml
camera_direction = camera_direction * 0.2;
quiver3(camera_pos(1), camera_pos(2), camera_pos(3), camera_direction(1), camera_direction(2), camera_direction(3), 0, 'LineWidth', 2);
%plot i_camera_direction
i_camera_direction = i_camera_direction * 0.2;
quiver3(camera_pos(1), camera_pos(2), camera_pos(3), i_camera_direction(1), i_camera_direction(2), i_camera_direction(3), 0, 'LineWidth', 2);
%plot j_camera_direction
j_camera_direction = j_camera_direction * 0.2;
quiver3(camera_pos(1), camera_pos(2), camera_pos(3), j_camera_direction(1), j_camera_direction(2), j_camera_direction(3), 0, 'LineWidth', 2);
% plot a camera at the origin
cam = plotCamera('Location', camera_pos(1:3), 'Orientation', cameraOrientation, 'Size', 0.1, 'Color', 'b');


function IAC = get_IAC(l_infs, vps, vp1s, vp2s, H)
%GET_IAC Returns the Image of the absolute conitc
%   l_infs set of imaged line at inf
%   vps set of vanishing points to be used with l_infs
%   vp1s set of vanishing points orthogonal to the point in the
%   corresponding pos in vp2
%   H is the homography in order to estimate the position of circular
%   points
%   IAC is the image of the absolute conic
% Assume w to have this form [a 0 b
%                             0 1 c
%                             b c d]

% matrix parametrization
syms a b c d;
omega = [a 0 b; 0 1 c; b c d];

X = []; % should be nxm matrix (n is ls size 2, m is 4)
Y = []; % should be n matrix of target values

% first add constraints on l_infs and vps
% 2 constraints for each couple
% [l_inf]x W vp = 0
eqn = [];
for ii = 1:size(l_infs,2)
    
    % first compute the element of l
    li = l_infs(:,ii);
    l1 = li(1,1);
    l2 = li(2,1);
    l3 = li(3,1);
    
    % vector product matrix
    lx = [0 -l3 l2; l3 0 -l1; -l2 l1 0];
    
    % get vp
    xi = vps(:,ii);
    
    eqn = [lx(1,:)*omega*xi == 0, lx(2,:)*omega*xi == 0];
    
end

% cast equations into matrix form
[A,y] = equationsToMatrix(eqn,[a,b,c,d]);
% concatenate matrices
X = [X;double(A)];
Y = [Y;double(y)];

% eqn contains all the equations
eqn = [];
% add constraints on vanishing points
for ii = 1:size(vp1s,2)
    % first compute the element of x
    vi = vp1s(:,ii);
    ui = vp2s(:,ii);
    
    % vp1' W vp2 = 0
    eqn = [eqn, vi.' * omega * ui == 0];
    
end

if size(eqn,2)>0
    % cast equations into matrix form
    [A,y] = equationsToMatrix(eqn,[a,b,c,d]);
    % concatenate matrices
    X = [X;double(A)];
    Y = [Y;double(y)];
end

% add constraints on homography
if size(H)>0
    % scaling factor needed in order to get an homogeneous matrix
    % get columns
    h1 = H(:,1);
    h2 = H(:,2);
    
    % first constraint: h1' w h2 = 0
    eq1 = h1.' * omega * h2 == 0;
    % second equation h1'wh1 = h2' w h2
    eq2 = h1.' * omega * h1 == h2.' * omega * h2;
    
    [A,y] = equationsToMatrix([eq1,eq2],[a,b,c,d]);
    A = double(A);
    y = double(y);
    
    % concatenate matrices
    X = [X;A];
    Y = [Y;y];
end

% fit a linear model without intercept
lm = fitlm(X,Y, 'y ~ x1 + x2 + x3 + x4 - 1');
% get the coefficients
W = lm.Coefficients.Estimate;

%W = X.'*X \ (X.'*Y)
% image of absolute conic
IAC = double([W(1,1) 0 W(2,1); 0 1 W(3,1); W(2,1) W(3,1) W(4,1)]);

end