% test_import_dxf
clc
clear
close all

pkg load geometry
pkg load miscellaneous

modelerPath = './Modeler_2D';
mesherPath = '../mesh2d';
meshWrapperPath = './Mesh_2D';
solverPath = './Solver_2D';
miscPath = './misc';

%%

dxfFileName = './sample_dxf.dxf';

% Global mesh properties

fMin = 1;
fMax = 3;
% 1. Maximum allowed mesh size with respect to wavelength.
% Take note that this takes into consideration the material
% properties.
relWLmeshMax = 0.33;

% 2. Bounding box additimesh.vert = on (Currently this is mesh
% properties, in the future this needs to belong to B.C.)
% This is given in relaitve terms of wavelength.
boundingBoxAddSpace = 0.125;

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

addpath(modelerPath);
addpath(genpath(mesherPath));
addpath(meshWrapperPath);
addpath(solverPath);
addpath(miscPath);

% Current DXF imported in mm
[header,tables,entities] = mod2D_importDXF(dxfFileName);

% Convert the read DXF file to mod2D "objects"
[dxfPolygons, ignoredEnts] = mod2D_convertEntityListToPolygons(entities);
fprintf(1,'Total %d entities ignored\n',numel(ignoredEnts));

% Reduce resolution to 50um
dxfPolygons = mod2D_reducePolygonResolution(dxfPolygons,200e-3);

% Subtract overlaps. This is to replace a manual macro.
dxfPolygons = mod2D_subtractOvelapingPolygons(dxfPolygons);

% Define mesh properties
meshProps = mesh2D_createMeshPropsStruct(...
              'relWLmeshMax',relWLmeshMax, ...
              'boundingBoxAddSpace',boundingBoxAddSpace, ...
              'algorithmType','delaunay',...
              'maxRadiusEdgeRatio',2.5);

% Define simulation properties
simProps = cem2D_createSimPropsStruct('fMin',fMin,'fMax',fMax,...
                                      'polarizationType','TE');

% Create the default material
defaultMaterial = cem2D_createMaterialDefs;
% Create the material list with only the default material
materialList = cem2D_addMaterialToList(defaultMaterial);

% Add a dielectric
material8p0 = cem2D_createMaterialDefs('name','plastic800','er',8);
materialList = cem2D_addMaterialToList(material8p0,materialList);

% Asign material properties to PCB circuit

materialAssignment = cell([1 numel(dxfPolygons)]);
materialAssignment(1:end) = {'plastic800'};

% Get default material properties
defaultMaterial = cem2D_getMaterialPropsFromName('default',materialList);

% Create bounding box (background)
bBox = mod2D_createBoundingBox(...
        dxfPolygons,...
        meshProps,...
        simProps,...
        defaultMaterial);

% Subtract all polygons from bounding box
for polIdx = 1:numel(dxfPolygons)
    bBox = mod2D_booleanOperation(bBox,dxfPolygons{polIdx},'subtract');
end

% Add bounding box to part list
polList = [{bBox}  dxfPolygons];
% Assign default material to bounding box
materialAssignment = [{'default'} materialAssignment];

% Draw some polygons yo!

figHdl = figure;
axHdl = axes;
mod2D_showPolygon(axHdl,polList{1},[1 1 1],[0 0 0]);
for polIdx = 1:numel(dxfPolygons)
    cPol = dxfPolygons{polIdx};
    mod2D_showPolygon(axHdl,cPol,[0.2 0.8 0.1],[0 0 0]);
end

warning off;

%%%%%%%%%
% Mesh! %
%%%%%%%%%

% 1. First phase of meas
initMesh = mesh2D_generateInitialMesh(...
            polList,...             % Entire polygon list. This is mainly to calculate the maximum\minimum edge length
            {},...                  % List of lines. Used for porst, mainly.
            materialAssignment,...  % Material assignments for initial LFS assignment
            materialList,...        % Corresponding material properties
            meshProps,...           % Darrens Mesher properties
            simProps);              % Simulation properties. This is to assign spatial meshing rules, most


% Plot initial mesh
figure('position',[240    165   1053    596]);
subplot(1,2,1);
% Draw bounding box
patch('faces',initMesh.tria(initMesh.tnum == 1,1:3),'vertices',initMesh.vert, ...
    'facecolor',[1,1,1], ...
    'edgecolor',[0,0,0]) ;
hold on;
axis image off;
% Plot the rest of the polygons
copperColor = [0.953 0.968 0.478];
for polIdx = 2:numel(polList)
    patch('faces',initMesh.tria(initMesh.tnum == polIdx,1:3),'vertices',initMesh.vert, ...
          'facecolor',copperColor, ...
          'edgecolor',[0,0,0]) ;
end

title(['MESH.: KIND=DELFRONT, |TRIA|=', ...
    num2str(size(initMesh.tria,1))]) ;

hold off;

% 2. Smooth the mesh
smoothMesh = mesh2D_smoothMesh(initMesh,meshProps);

subplot(1,2,2);
patch('faces',smoothMesh.tria(smoothMesh.tnum == 1,1:3),'vertices',smoothMesh.vert, ...
    'facecolor',[1.,1.,1.], ...
    'edgecolor',[0,0,0]) ;
hold on;
axis image off;
for polIdx = 2:numel(polList)
    patch('faces',initMesh.tria(initMesh.tnum == polIdx,1:3),'vertices',initMesh.vert, ...
          'facecolor',copperColor, ...
          'edgecolor',[0,0,0]) ;
end

title(['MESH-Smoothed.: KIND=DELFRONT, |TRIA|=', ...
    num2str(size(smoothMesh.tria,1))]) ;

hold off;

warning on;

rmpath(modelerPath);
rmpath(genpath(mesherPath));
rmpath(meshWrapperPath);
rmpath(solverPath);
rmpath(miscPath);
