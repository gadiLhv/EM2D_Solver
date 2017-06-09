% test_geometry_package.m
clc
clear
close all

pol1 = [2 2; 6 2; 6 6; 2 6; 2 2; NaN NaN; 3 3; 3 5; 5 5; 5 3; 3 3];
pol1a = [2 6; 2 2; 3 3; 3 5; 5 5; 5 3; 3 3; 2 2; 6 2; 6 6; 2 6];
pol2 = [1 2; 7 4; 4 7; 1 2; NaN NaN; 2.5 3; 5.5 4; 4 5.5; 2.5 3];
pol2a = [1 2; 7 4; 4 7; 1 2; 2.5 3; 4 5.5; 5.5 4; 2.5 3; 1 2];
lw = 2;

subplot (2, 6, [2 3])
patch (pol1a(:, 1), pol1a(:, 2), 'facecolor', 'c', 'edgecolor', 'k', 'linewidth', lw);
axis image
grid on
title ("1. Subject polygon")

subplot (2, 6, [4 5])
patch (pol1a(:, 1), pol1a(:, 2), 'facecolor', 'c', 'edgecolor', 'none');
hold on
patch (pol2a(:, 1), pol2a(:, 2), 'facecolor', 'y', 'edgecolor', 'b', 'linewidth', lw);
axis image
grid on
title "2. Clip polygon"

op   = {"Sub -clip", "AND / Intersection", "Exclusive OR", "OR / Union"};

for i=1:numel(op)
  subplot (6, 4, [12 16]+i);
  [opol, npol] = clipPolygon(pol1, pol2, i-1);
  opol = polygon2patch (opol);
  patch (pol1a(:, 1), pol1a(:, 2), 'facecolor', 'c', 'edgecolor', 'none');
  hold on
  patch (pol2a(:, 1), pol2a(:, 2), 'facecolor', 'y', 'edgecolor', 'none');
  patch (opol(:,1),opol(:,2), 'facecolor', 'g', 'edgecolor', 'r', ...
       'linewidth', lw, 'erasemode', 'xor');
  axis image
  grid on
  title (sprintf("%d. %s", i+2, op{i}));
  axis off
end

subplot (10, 4, 37);
[opol, npol] = clipPolygon(pol2, pol1, 0);
opol = polygon2patch (opol);
patch (pol1a(:, 1), pol1a(:, 2), 'facecolor', 'c', 'edgecolor', 'none');
hold on
patch (pol2a(:, 1), pol2a(:, 2), 'facecolor', 'y', 'edgecolor', 'none');
patch (opol(:,1),opol(:,2), 'facecolor', 'g', 'edgecolor', 'r', ...
     'linewidth', lw, 'erasemode', 'xor');
axis image
grid on
axis off
title "7. Clip - sub";