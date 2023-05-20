% test_interpolant_functions.m
clc
clear
close all

warning off;

r1 = [0 0];
r2 = [1 0];
r3 = [0.7 0.75];


%%

vert_m = [r1 ; r2 ; r3];
[a,b,c,Det] = cem2D_createInterpolantCoeffs1st(vert_m,[1 2 3]);
l1 = sqrt(sum((r2 - r1).^2,2));
l2 = sqrt(sum((r3 - r2).^2,2));
l3 = sqrt(sum((r1 - r3).^2,2));

Li = @(i,x,y) (1/(2*Det))*(a(i) + b(i)*x + c(i)*y);
dLidx = @(i,x,y) b(i);
dLidy = @(i,x,y) c(i);

N1x = @(x,y) l1*(Li(1,x,y)*b(2) - Li(2,x,y)*b(1));
N1y = @(x,y) l1*(Li(1,x,y)*c(2) - Li(2,x,y)*c(1));

N2x = @(x,y) l2*(Li(2,x,y)*b(3) - Li(3,x,y)*b(2));
N2y = @(x,y) l2*(Li(2,x,y)*c(3) - Li(3,x,y)*c(2));

N3x = @(x,y) l3*(Li(3,x,y)*b(1) - Li(1,x,y)*b(3));
N3y = @(x,y) l3*(Li(3,x,y)*c(1) - Li(1,x,y)*c(3));




%%

% Create grid

nx = 11;
ny = 11;

[xm,ym] = meshgrid(...
    linspace(min(vert_m(:,1)),max(vert_m(:,1)),nx),...
    linspace(min(vert_m(:,2)),max(vert_m(:,2)),ny));

bIn = inpolygon(xm(:), ym(:),vert_m(:,1), vert_m(:,2));

figure('position', [300    200   1202    770]);
subplot(2,2,1);
polHdl = patch (vert_m(:,1),vert_m(:,2), ...
              'facecolor', [1 1 1], ...
              'edgecolor', [0 0 0]);




hold on;
quiver(xm(bIn),ym(bIn),N1x(xm(bIn),ym(bIn)),N1y(xm(bIn),ym(bIn)));
hold off;

subplot(2,2,2);
polHdl = patch (vert_m(:,1),vert_m(:,2), ...
              'facecolor', [1 1 1], ...
              'edgecolor', [0 0 0]);



hold on;
quiver(xm(bIn),ym(bIn),N2x(xm(bIn),ym(bIn)),N2y(xm(bIn),ym(bIn)));
hold off;

subplot(2,2,3);
polHdl = patch (vert_m(:,1),vert_m(:,2), ...
              'facecolor', [1 1 1], ...
              'edgecolor', [0 0 0]);



hold on;
quiver(xm(bIn),ym(bIn),N3x(xm(bIn),ym(bIn)),N3y(xm(bIn),ym(bIn)));
hold off;

warning on;
