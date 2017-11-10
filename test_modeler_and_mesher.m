% test_modeler_and_mesher
clc
clear
close all

modelerPath = './Modeler_2D';
mesherPath = '/home/gadi/Repositories/mesh2d';
meshWrapperPath = './Mesh_2D';
solverPath = './Solver_2D';


%% Mesher parameters
% Set printout/no printout
outputMesh = true(1);

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

% Create default material


addpath(modelerPath);
addpath(mesherPath);
addpath(meshWrapperPath);
addpath(solverPath);

% Define mesh properties
meshProps = mesh2D_createMeshPropsStruct(...
              'relWLmeshMax',relWLmeshMax, ...
              'boundingBoxAddSpace',boundingBoxAddSpace);

% Define simulation properties
simProps = cem2D_createSimPropsStruct('fMin',fMin,'fMax',fMax,...
                                      'polarizationType','TE');
              
% Create the default material
defaultMaterial = cem2D_createMaterialDefs;
% Create the material list with only the default material
materialList = cem2D_addMaterialToList(defaultMaterial);

% Add three materials
material3p1 = cem2D_createMaterialDefs('name','plastic310','er',3.1);
material5p5 = cem2D_createMaterialDefs('name','plastic550','er',5.5);
material7p2 = cem2D_createMaterialDefs('name','plastic720','er',7.2);

materialList = cem2D_addMaterialToList(material3p1,materialList);
materialList = cem2D_addMaterialToList(material5p5,materialList);
materialList = cem2D_addMaterialToList(material7p2,materialList);

%% Create the geometry


% Geometry 1
pol1 = mod2D_createRectangleStruct([-2 -2]*10,[2 -1]*10);
pol2 = mod2D_createRectangleStruct([-2 -1]*10,[2 0]*10);
pol3 = mod2D_createRectangleStruct([-2 0]*10,[2 1]*10);
rect1 = mod2D_createRectangleStruct([-1 -1.5]*10,[1 0.5]*10);
pol1 =  mod2D_booleanOperation(pol1,rect1,'subtract');
pol2 =  mod2D_booleanOperation(pol2,rect1,'subtract');
pol3 =  mod2D_booleanOperation(pol3,rect1,'subtract');

% Geometery 2
%rect1 = mod2D_createRectangleStruct([-3 -3],[3 3]);
%rect2 = mod2D_createRectangleStruct([-2.5 -2.5],[2.5 2.5]);
%pol1 = mod2D_booleanOperation(rect1,rect2,'subtract');
%rect1 = mod2D_createRectangleStruct([-2 -2],[2 2]);
%rect2 = mod2D_createRectangleStruct([-1 -1],[1 1]);
%pol2 = mod2D+_booleanOperation(rect1,rect2,'subtract');
%pol1 = mod2D_booleanOperation(pol1,pol2,'add');
%rect1 = mod2D_createRectangleStruct([2 -0.5],[2.5 1.5]);
%pol1 = mod2D_booleanOperation(pol1,rect1,'add');
%pol2 = mod2D_createRectangleStruct([-2.5 1],[-2 2]); % Bad at handling duplicate nodes. Remove them while creating connectivty list
%pol3 = mod2D_createPolygonStruct([-0.7 0.4 1],[-0.6 -0.8 1]);

% Create polygon list
polList = {pol1,pol2,pol3};

% Assign materials,er,mr
materialAssignment = {'plastic310','plastic550','plastic720'};

% Get default material properties
defaultMaterial = cem2D_getMaterialPropsFromName('default',materialList);

% Create bounding box (background)
bBox = mod2D_createBoundingBox(...
        polList,...
        meshProps,...
        simProps,...
        defaultMaterial);

% Assign materials to polygons


% Add all polygons before subtracting from the background
if(numel(polList) < 2)
  allPols = polList{1};
elseif(numel(polList) >= 2)
  allPols = mod2D_booleanOperation(polList{1},polList{2},'add');
  if(numel(polList) > 2)
    for polIdx = 3:numel(polList)
      allPols = mod2D_booleanOperation(allPols,polList{polIdx},'add');
    end
  end
end

% Subtract all polygons from background
bBox = mod2D_booleanOperation(bBox,allPols,'subtract');

% Add bounding box to part list
polList = [{bBox}  polList];

% Assign default material to bounding box
materialAssignment = [{'default'} materialAssignment];

[node,edge,face] = mod2D_polygonToFaceList(polList,meshProps.stitchingTolerance);


%%%%%%%%%
% Mesh! %
%%%%%%%%%

% 1. First phase of meas
initMesh = mesh2D_generateInitialMesh(...
            polList,...             % Entire polygon list. This is mainly to calculate the maximum\minimum edge length
            materialAssignment,...  % Material assignments for initial LFS assignment
            materialList,...        % Corresponding material properties
            meshProps,...           % Darrens Mesher properties
            simProps);              % Simulation properties. This is to assign spatial meshing rules, most

% 2. Smooth the mesh
smoothMesh = mesh2D_smoothMesh(initMesh);
            
% Plot initial mesh
figure('position',[240    165   1053    596]);
subplot(1,2,1);
patch('faces',initMesh.tria(initMesh.tnum == 1,1:3),'vertices',initMesh.vert, ...
    'facecolor',[1.,1.,1.], ...
    'edgecolor',[0,0,0]) ;
hold on; 
axis image off;
randColor = rand([1 3]);
patch('faces',initMesh.tria(initMesh.tnum == 2,1:3),'vertices',initMesh.vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
randColor = rand([1 3]);
patch('faces',initMesh.tria(initMesh.tnum == 3,1:3),'vertices',initMesh.vert, ...
    'facecolor',randColor,...
    'edgecolor',[0,0,0]) ;
randColor = rand([1 3]);
patch('faces',initMesh.tria(initMesh.tnum == 4,1:3),'vertices',initMesh.vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
    
edgeVect1 = initMesh.vert(initMesh.etri(:,1),:);
edgeVect2 = initMesh.vert(initMesh.etri(:,2),:);
edgeVectX = [edgeVect1(:,1) edgeVect2(:,1)].';
edgeVectY = [edgeVect1(:,2) edgeVect2(:,2)].';
plot(edgeVectX,edgeVectY,'-','linewidth',3);
    
title(['MESH.: KIND=DELFRONT, |TRIA|=', ...
    num2str(size(initMesh.tria,1))]) ;
    
hold off;

subplot(1,2,2);
patch('faces',smoothMesh.tria(smoothMesh.tnum == 1,1:3),'vertices',smoothMesh.vert, ...
    'facecolor',[1.,1.,1.], ...
    'edgecolor',[0,0,0]) ;
hold on; 
axis image off;
randColor = rand([1 3]);
patch('faces',smoothMesh.tria(smoothMesh.tnum == 2,1:3),'vertices',smoothMesh.vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
randColor = rand([1 3]);
patch('faces',smoothMesh.tria(smoothMesh.tnum == 3,1:3),'vertices',smoothMesh.vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
randColor = rand([1 3]);
patch('faces',smoothMesh.tria(smoothMesh.tnum == 4,1:3),'vertices',smoothMesh.vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
title(['MESH-Smoothed.: KIND=DELFRONT, |TRIA|=', ...
    num2str(size(smoothMesh.tria,1))]) ;
hold off;

rmpath(modelerPath);
rmpath(mesherPath);
rmpath(meshWrapperPath);
rmpath(solverPath);