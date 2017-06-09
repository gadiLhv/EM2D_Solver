% test_polygons_and_squares
clc
clear
close all

modelerPath = './Modeler_2D';
%%

addpath(modelerPath);

rect1 = mod2D_createRectangleStruct([0.5 0.5],[2.5 2.5]);
rect2 = mod2D_createRectangleStruct([1.5 1.5],[2 2]);

holeRect = mod2D_booleanOperation(rect1,rect2,'subtract');

figHdl = figure;
axHdl = axes;
axis([-4 4 -4 4]);

mod2D_showPolygon(axHdl,holeRect,[0.2 0.8 0.1],[0 0 0]);

tri = mod2D_createPolygonStruct([1.7 2 0],[1.7 -1 0]);
newPoly = mod2D_booleanOperation(holeRect,tri,'add');
pause(1);

clf(figHdl);
axHdl = axes;
axis([-4 4 -4 4]);
mod2D_showPolygon(axHdl,newPoly,[0.2 0.8 0.1],[0 0 0]);

rmpath(modelerPath);