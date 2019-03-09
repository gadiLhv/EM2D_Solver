% test_import_dxf
clc
clear
close all

modelerPath = './Modeler_2D';

addpath(modelerPath);

[header,tables,entities] = mod2D_importDXF('./sample_dxf.dxf');

[mod2Dpolygons, ignoredEnts] = mod2D_convertEntityListToPolygons(entities);
fprintf(1,'Total %d entities ignored\n',numel(ignoredEnts));

mod2Dpolygons = mod2D_reducePolygonResolution(mod2Dpolygons,50e-3);

mod2Dpolygons = mod2D_subtractOvelapingPolygons(mod2Dpolygons);


% Draw some polygons yo!

figHdl = figure;
axHdl = axes;
for polIdx = 1:numel(mod2Dpolygons)
    cPol = mod2Dpolygons{polIdx};
    mod2D_showPolygon(axHdl,cPol,[0.2 0.8 0.1],[0 0 0]);
end

rmpath(modelerPath);