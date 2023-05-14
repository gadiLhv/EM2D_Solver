% test_mesh.m

clc
clear
close all

pkg load geometry

modelerPath = './Modeler_2D';
mesherPath = '../mesh2d';
meshWrapperPath = './Mesh_2D';
solver2DPath = './Solver_2D';
solver1DPath = './Solver_1D';
miscPath = './misc';


dxfFile = 'C:\Users\Gadi_Lahav\Documents\KiCAD\RF_Test_Board\RF_Test_Board\Multi_Layer\Multi_Layer-F_Cu.dxf';

%%

pkg('load','miscellaneous');
addpath(modelerPath);
addpath(genpath(mesherPath));
addpath(meshWrapperPath);
addpath(solver2DPath);
addpath(solver1DPath);
addpath(miscPath);


[~,~,dxfEnts] = mod2D_importDXF(dxfFile);

[pols,ignoredEnts] = mod2D_convertEntityListToPolygons(dxfEnts);

pol = pols{1};

cx = mean(pol.x);
cy = mean(pol.y);

pol.x = pol.x - cx;
pol.y = pol.y - cy;
bBox = max(max(abs(pol.x)),max(abs(pol.y)))*[-1 1 -1 1];

polList = {pol};
materialAssignment = {'default'};

% Create the materials for this simulation
material_default = cem2D_createMaterialDefs;
% Create the material list with only the default material
materialList = cem2D_addMaterialToList(material_default);

% Simulation properties
simProps = cem2D_createSimPropsStruct(...
                'fMin',0.5,...
                'fMax',1,...
                'polarizationType','TE');


% Global mesh properties
meshProps = mesh2D_createMeshPropsStruct(...
              'relWLmeshMax',0.33, ...
              'boundingBoxAddSpace',1, ...
              'algorithmType','delfront',...
              'maxRadiusEdgeRatio',1.5);



meshData = mesh2D_generateInitialMesh(...
            polList,...             % Entire polygon list. This is mainly to calculate the maximum\minimum edge length
            {},...                  % List of lines. Used for porst, mainly.
            materialAssignment,...  % Material assignments for initial LFS assignment
            materialList,...        % Corresponding material properties
            meshProps,...           % Darrens Mesher properties
            simProps);              % Simulation properties. This is to assign spatial meshing rules, most



meshData = mesh2D_smoothMesh(meshData,meshProps);

figure('position',[300   202   722   632]);

%patch('faces',meshData.tria(meshData.tnum == 1,1:3),'vertices',meshData.vert, ...
%    'facecolor',[234,207,181]/255, ...
%    'edgecolor',[0,0,0]) ;

patch('faces',meshData.tria(meshData.tnum == 1,1:3),'vertices',meshData.vert, ...
    'facecolor',[1 1 1], ...
    'edgecolor',[0,0,0]) ;

axis(bBox);
axis('square');
axis off;

pause(1);
print('vivaldi_mesh.svg','-dsvg');



