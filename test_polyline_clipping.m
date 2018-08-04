% test_polyline_clipping.m
clc
clear
close all


%sline = [0, 6.5; 1.25, 4; 1.25, 0; NaN, NaN; 0.25, 7; 1.75, 4; 1.75, 0];
%for ii=1:10
%  sline = [sline; [NaN NaN]; [ii/2+0.25, 7; ii/2+1.75, 4; ii/2+1.75, 0]];
%endfor
sline = [5/2+0.25, 7; 5/2+1.75, 4; 5/2+1.75, 0];
pol2a = [1 2; 7 4; 4 7; 1 2; 2.5 3; 4 5.5; 5.5 4; 2.5 3; 1 2];
figure;
subplot(2,1,1);
plot(sline(:, 1), sline(:, 2));
hold on;
patch(pol2a(:, 1), pol2a(:, 2), 'facecolor', 'y', 'edgecolor', 'b', 'linewidth', 2);
hold off;
grid on;
axis equal;
title('Original polygon and lines');

subplot(2,2,3);
% "AND" operation
patch(pol2a(:, 1), pol2a(:, 2), 'facecolor', 'y', 'edgecolor', 'b', 'linewidth', 2);
hold on;
[olin, nlin] = clipPolyline (sline, pol2a, 1);
plot (olin(:, 1), olin(:, 2), 'r', 'linewidth', 3);
hold off;
grid on;
axis equal;
title('"AND"\Intersect');

subplot(2,2,4);
% Subtract operation
patch(pol2a(:, 1), pol2a(:, 2), 'facecolor', 'y', 'edgecolor', 'b', 'linewidth', 2);
hold on;
[olin, nlin] = clipPolyline(sline, pol2a, 0);
plot (olin(:, 1), olin(:, 2), 'g', 'linewidth', 3);
hold off;
grid on;
axis equal;
title ('Clip\Subtract');

