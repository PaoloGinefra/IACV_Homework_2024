close all;
clc;
clear;

import java.util.LinkedList;

%% Load Image
im = imread('../../Assignment/Homework Image.jpg');

figure(1);
imshow(im);
alpha(0.5);
hold on;
drawnow;

[x, y] = ginput(1);
x = round(x);
y = round(y);
disp(['Selected pixel: ', num2str(x), ', ', num2str(y)]);

threshold = 5;
frontier = LinkedList();
frontier.add([x, y]);

visited = zeros(size(im, 1), size(im, 2));


while ~frontier.isEmpty()
    current = frontier.remove();
    x = current(1);
    y = current(2);
    
    if visited(y, x) == 1
        continue;
    end
    
    visited(y, x) = 1;
    
    figure(1);
    plot(x, y, 'r.');
    drawnow;
    
    deltas = [-1, 0; 1, 0; 0, -1; 0, 1];
    
    for i = 1:size(deltas, 1)
        x1 = x + deltas(i, 1);
        y1 = y + deltas(i, 2);
        
        if x1 < 1 || x1 > size(im, 2) || y1 < 1 || y1 > size(im, 1)
            continue;
        end
        
        if visited(y1, x1) == 1
            continue;
        end
        
        if ComputeDistance(im, x, y, x1, y1) < threshold
            frontier.add([x1, y1]);
        end
    end
    
end

hold off;



function [computedDist] = ComputeDistance(im, x, y, x1, y1)
p1 = im(y, x, :);
p2 = im(y1, x1, :);
computedDist = sqrt(sum((p1 - p2).^2));
end