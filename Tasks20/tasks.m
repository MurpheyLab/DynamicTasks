% tasks.m
% script to generate tasks for ball-in-bowl experiment
% last edited 09.05.19 by Kyra Rudy
close all

%% CLUSTERS
clear
filename = 'task.csv';

% set parameters for cluster count and size
clusters = 5;
cluster_size = 4;
spread = 0.05;

dots = clusters * cluster_size;
goals = zeros(dots, 3);

% randomly populate clusters of points
for i = 1:clusters
    goals(1 + (i-1)*cluster_size, 1:2) = [rand(1), rand(1)];
    for j = 1:cluster_size-1
        goals(1 + (i-1)*cluster_size + j, 1:2) = goals(1 + (i-1)*cluster_size, 1:2) + ...
                                                 spread .* [randn(1), randn(1)];
    end
end

% make sure none of the goals are out of range
for i = 1:dots
    for j = 1:2
        if goals(i,j) > 1
            goals(i,j) = 0.99;
        elseif goals(i,j) < 0
            goals(i,j) = 0.01;
        end
    end
end

figure()
plot(goals(:,1), goals(:,2), '*');

% save to csv file
writematrix(goals, filename);

%% BOWTIE
clear

dots = 20;
goals = zeros(dots, 3);

% horizontal line
goals(1:6,1) = linspace(0.15, 0.85, 6)';
goals(1:6,2) = 0.5 .* ones(6,1);

% vertical line
goals(7:12,2) = linspace(0.15, 0.85, 6)';
goals(7:12,1) = 0.5 .* ones(6,1);

% diagonals
goals(13:16,1) = linspace(0.5, 0.85, 4)';
goals(13:16,2) = linspace(0.85, 0.5, 4)';
goals(17:dots,1) = linspace(0.15, 0.5, 4)';
goals(17:dots,2) = linspace(0.5, 0.15, 4)';

figure()
plot(goals(:,1), goals(:,2), '*');

% save to csv file
writematrix(goals, 'bowtie_task.csv');

%% SPIRAL
clear

% spiral parameters
pos = [0.5 0.5;    % startpoint
       0.5 0.9];  % endpoint
% pos = [0.4 0.4;    % startpoint
%        0.7 0.7];  % endpoint
nturns = 2;    % number of turns (integer value)
dots = 20;
goals = zeros(dots, 3);

% define spiral
dp = diff(pos,1,1);
R = hypot(dp(1), dp(2));
phi0 = atan2(dp(2), dp(1));
phi = linspace(0, nturns*2*pi, dots);
r = linspace(0, R, numel(phi));

x = pos(1,1) + r .* cos(phi + phi0);
y = pos(1,2) + r  .* sin(phi + phi0);

goals(1:dots,1) = x';
goals(1:dots,2) = y';

figure()
plot(goals(:,1), goals(:,2), '*'); 

% save to csv file
writematrix(goals, 'spiral_task.csv');

%% EXAMPLE TASK (COMBO)
clear
filename = 'ex_task.csv';

% set parameters for cluster count and size
clusters = 2;
cluster_size = 5;
spread = 0.05;

dots = 2 * clusters * cluster_size;
goals = zeros(dots, 3);

% randomly populate clusters of points
for i = 1:clusters
    goals(1 + (i-1)*cluster_size, 1:2) = [0.5*rand(1), rand(1)];
    for j = 1:cluster_size-1
        goals(1 + (i-1)*cluster_size + j, 1:2) = goals(1 + (i-1)*cluster_size, 1:2) + ...
                                                 spread .* [randn(1), randn(1)];
    end
end

% make sure none of the goals are out of range
for i = 1:dots/2
    for j = 1:2
        if goals(i,j) > 1
            goals(i,j) = 0.99;
        elseif goals(i,j) < 0
            goals(i,j) = 0.01;
        end
    end
end

line = linspace(0.2,0.8,10);
goals((dots/2 + 1):dots,2) = line';
for i = (dots/2 + 1):dots
   goals(i,1) = [0.75]; 
end

figure()
plot(goals(:,1), goals(:,2), '*');

% save to csv file
writematrix(goals, filename);



