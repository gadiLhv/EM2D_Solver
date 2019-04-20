% test_modeler_and_mesher
clc
clear
close all

modelerPath = './Modeler_2D';
mesherPath = '/home/gadi/Repositories/mesh2d';
meshWrapperPath = './Mesh_2D';
solverPath = './Solver_2D';
miscPath = './misc';

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
addpath(genpath(mesherPath));
addpath(meshWrapperPath);
addpath(solverPath);
addpath(miscPath);

% Define mesh properties
meshProps = mesh2D_createMeshPropsStruct(...
              'relWLmeshMax',relWLmeshMax, ...
              'boundingBoxAddSpace',boundingBoxAddSpace, ...
              'algorithmType','delaunay',... delaunay delfront
              'maxRadiusEdgeRatio',1.25);

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
rect1 = mod2D_createRectangleStruct([-1 -0.5]*10,[1 0.5]*10);
rect2 = mod2D_createRectangleStruct([-0.5 -1.75]*10,[0.5 -1.25]*10);
pol1 =  mod2D_booleanOperation(pol1,rect2,'subtract');
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

% Create edge port line
line1 = mod2D_createLineStruct([-1.5 -1.5]*10,[0.3 0.8]*10);

% Create polygon list
polList = {pol1,pol2,pol3};
lineList = {line1};

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

%[node,edge,face] = mod2D_polygonToFaceList(polList,meshProps.stitchingTolerance);

%%%%%%%%%
% Mesh! %
%%%%%%%%%

% 1. First phase of meas
initMesh = mesh2D_generateInitialMesh(...
            polList,...             % Entire polygon list. This is mainly to calculate the maximum\minimum edge length
            lineList,...            % List of lines. Used for porst, mainly.
            materialAssignment,...  % Material assignments for initial LFS assignment
            materialList,...        % Corresponding material properties
            meshProps,...           % Darrens Mesher properties
            simProps);              % Simulation properties. This is to assign spatial meshing rules, most

% 2. Smooth the mesh
finalMesh = mesh2D_smoothMesh(initMesh,meshProps);
            
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
plot(edgeVectX,edgeVectY,'-','linewidth',2,'color',[0 0 0]);

edgeVect1 = initMesh.vert(initMesh.ltri(:,1),:);
edgeVect2 = initMesh.vert(initMesh.ltri(:,2),:);
edgeVectX = [edgeVect1(:,1) edgeVect2(:,1)].';
edgeVectY = [edgeVect1(:,2) edgeVect2(:,2)].';
%plot(edgeVectX,edgeVectY,'-','linewidth',2,'color',[1 0 0]);
for segIdx = 1:size(edgeVectX,2)
  plot([edgeVectX(1,segIdx) edgeVectX(2,segIdx)],...
       [edgeVectY(1,segIdx) edgeVectY(2,segIdx)],...
       '-','linewidth',2,'color',[1 0 0]);  
end
    
title(['MESH.: KIND=DELFRONT, |TRIA|=', ...
    num2str(size(initMesh.tria,1))]) ;
    
hold off;

subplot(1,2,2);
patch('faces',finalMesh.tria(finalMesh.tnum == 1,1:3),'vertices',finalMesh.vert, ...
    'facecolor',[1.,1.,1.], ...
    'edgecolor',[0,0,0]) ;
hold on; 
axis image off;
randColor = rand([1 3]);
patch('faces',finalMesh.tria(finalMesh.tnum == 2,1:3),'vertices',finalMesh.vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
randColor = rand([1 3]);
patch('faces',finalMesh.tria(finalMesh.tnum == 3,1:3),'vertices',finalMesh.vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
randColor = rand([1 3]);
patch('faces',finalMesh.tria(finalMesh.tnum == 4,1:3),'vertices',finalMesh.vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;

edgeVect1 = finalMesh.vert(finalMesh.etri(:,1),:);
edgeVect2 = finalMesh.vert(finalMesh.etri(:,2),:);
edgeVectX = [edgeVect1(:,1) edgeVect2(:,1)].';
edgeVectY = [edgeVect1(:,2) edgeVect2(:,2)].';
plot(edgeVectX,edgeVectY,'-','linewidth',2,'color',[0 0 0]);
    
edgeVect1 = finalMesh.vert(finalMesh.ltri(:,1),:);
edgeVect2 = finalMesh.vert(finalMesh.ltri(:,2),:);
edgeVectX = [edgeVect1(:,1) edgeVect2(:,1)].';
edgeVectY = [edgeVect1(:,2) edgeVect2(:,2)].';
%plot(edgeVectX,edgeVectY,'-','linewidth',2,'color',[1 0 0]);
for segIdx = 1:size(edgeVectX,2)
  plot([edgeVectX(1,segIdx) edgeVectX(2,segIdx)],...
       [edgeVectY(1,segIdx) edgeVectY(2,segIdx)],...
       '-','linewidth',2,'color',[1 0 0]);
end

title(['MESH-Smoothed.: KIND=DELFRONT, |TRIA|=', ...
    num2str(size(finalMesh.tria,1))]) ;
    
hold off;

rmpath(modelerPath);
rmpath(genpath(mesherPath));
rmpath(meshWrapperPath);
rmpath(solverPath);