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
triTriplets = [1 2 3];
edgeTriplets = [1 2 3];

%%

% Create grid

nx = 15;
ny = 15;

[xm,ym] = meshgrid(...
    linspace(min(vert_m(:,1)),max(vert_m(:,1)),nx),...
    linspace(min(vert_m(:,2)),max(vert_m(:,2)),ny));




figure('position', [300    200   1202    770]);
subplot(2,2,1);
polHdl = patch (vert_m(:,1),vert_m(:,2), ...
              'facecolor', [1 1 1], ...
              'edgecolor', [0 0 0]);


fI_1 = cem2D_vectorElementInterp(vert_m,triTriplets,edgeTriplets,[1 0 0].');
fI_2 = cem2D_vectorElementInterp(vert_m,triTriplets,edgeTriplets,[0 1 0].');
fI_3 = cem2D_vectorElementInterp(vert_m,triTriplets,edgeTriplets,[0 0 1].');


N1 = fI_1(xm,ym);
N2 = fI_2(xm,ym);
N3 = fI_3(xm,ym);

xm = xm(:);
ym = ym(:);

hold on;
quiver(xm,ym,N1(:,1),N1(:,2));
hold off;

subplot(2,2,2);
polHdl = patch (vert_m(:,1),vert_m(:,2), ...
              'facecolor', [1 1 1], ...
              'edgecolor', [0 0 0]);



hold on;
quiver(xm,ym,N2(:,1),N2(:,2));
hold off;

subplot(2,2,3);
polHdl = patch (vert_m(:,1),vert_m(:,2), ...
              'facecolor', [1 1 1], ...
              'edgecolor', [0 0 0]);



hold on;
quiver(xm,ym,N3(:,1),N3(:,2));
hold off;

warning on;
