% test_GCPW_portmodes.m

clc
clear
close all


modelerPath = './Modeler_2D';
mesherPath = '../mesh2d';
meshWrapperPath = './Mesh_2D';
solver2DPath = './Solver_2D';
solver1DPath = './Solver_1D';
miscPath = './misc';

warning('off')
pkg('load','matgeom');
pkg('load','geometry');
pkg('load','miscellaneous');
addpath(modelerPath);
addpath(genpath(mesherPath));
addpath(meshWrapperPath);
addpath(solver2DPath);
addpath(solver1DPath);
addpath(miscPath);

% Simulation properties
simProps = cem2D_createSimPropsStruct(...
                'fMin',2,...
                'fMax',3,...
                'polarizationType','TE');

% Global mesh properties
meshProps = mesh2D_createMeshPropsStruct(...
              'relWLmeshMax',0.33, ...
              'boundingBoxAddSpace',0, ...
              'algorithmType','delfront',...
              'maxRadiusEdgeRatio',1.25);

              
% GCPW parameters

gcpw_W_mm = 0.5;                        % Conductor width
gcpw_g_mm = 0.5;                        % gap between ground and conductor.
gcpw_t_mm = 18/1e3;                   % Conductor thickness
gcpw_sub_thick_mm = 200/1e3;          % Substrate thickness
port_W_mm = 10;                       % Width of the port
port_H_mm = 5;                        % Height above copper thickness

topGround = mod2D_createRectangleStruct([-0.5*port_W_mm gcpw_sub_thick_mm],[0.5*port_W_mm gcpw_sub_thick_mm+gcpw_t_mm]);
botGround = mod2D_createRectangleStruct([-0.5*port_W_mm -gcpw_t_mm],[0.5*port_W_mm 0]);
substrate = mod2D_createRectangleStruct([-0.5*port_W_mm 0],[0.5*port_W_mm gcpw_sub_thick_mm]);
airBox = mod2D_createRectangleStruct([-0.5*port_W_mm -gcpw_t_mm],[0.5*port_W_mm gcpw_sub_thick_mm+gcpw_t_mm+port_H_mm]);

gcpwGndCutout1 = mod2D_createRectangleStruct([-0.5*gcpw_W_mm-gcpw_g_mm gcpw_sub_thick_mm],[-0.5*gcpw_W_mm gcpw_sub_thick_mm+gcpw_t_mm]);
gcpwGndCutout2 = mod2D_createRectangleStruct([0.5*gcpw_W_mm gcpw_sub_thick_mm],[0.5*gcpw_W_mm+gcpw_g_mm gcpw_sub_thick_mm+gcpw_t_mm]);

topGround = mod2D_booleanOperation(topGround,gcpwGndCutout1,'subtract');
topGround = mod2D_booleanOperation(topGround,gcpwGndCutout2,'subtract');

%%%%%%%%%%%%%%%%%%%%
% Assign materials %
%%%%%%%%%%%%%%%%%%%%

% Create the materials for this simulation
material_default = cem2D_createMaterialDefs;
material_PEC = cem2D_createMaterialDefs('name','PEC','type','PEC');
material_substrate = cem2D_createMaterialDefs(...
    'name','FR-4',...
    'er',4.3,...
    'tand_e',0.02,...
    'f_m',2);

% Create the material list with only the default material
materialList = cem2D_addMaterialToList(material_default);
materialList = cem2D_addMaterialToList(material_PEC,materialList);
materialList = cem2D_addMaterialToList(material_substrate,materialList);

% Create polygon list 
polList = {substrate,topGround,botGround};
% Assign materials per-shape
materialAssignment = {'FR-4','PEC','PEC'};


% Subtract everything from airbox
airBox = mod2D_booleanOperation(airBox,topGround,'subtract');
airBox = mod2D_booleanOperation(airBox,substrate,'subtract');
airBox = mod2D_booleanOperation(airBox,botGround,'subtract');

polList = [{airBox} polList];
materialAssignment = [{'default'},materialAssignment];

figHdl = figure('position',[108    203   1174    417]);
axHdl = axes;
mod2D_showPolygon(axHdl,botGround,[204 204 0]/255,[0 0 0]);
mod2D_showPolygon(axHdl,substrate,[150 255 150]/255,[0 0 0]);
mod2D_showPolygon(axHdl,topGround,[204 204 0]/255,[0 0 0]);
mod2D_showPolygon(axHdl,airBox,[255 255 255]/255,[0 0 0]);
set(gca,'fontsize',14);
xlabel('x [mm]','fontsize',16);
ylabel('y [mm]','fontsize',16);



% Test mesh of bounding box
meshData = mesh2D_generateInitialMesh(...
            polList,...             % Entire polygon list. This is mainly to calculate the maximum\minimum edge length
            {},...                  % List of lines. Used for porst, mainly.
            materialAssignment,...  % Material assignments for initial LFS assignment
            materialList,...        % Corresponding material properties
            meshProps,...           % Darrens Mesher properties
            simProps);              % Simulation properties. This is to assign spatial meshing rules, most

% Smooth the mesh
meshData = mesh2D_smoothMesh(meshData,meshProps);

hold(axHdl,'on');
patch('faces',meshData.tria(meshData.tnum == 1,1:3),'vertices',meshData.vert, ...
    'facecolor','none', ...
    'edgecolor',[0,0,0]) ;

patch('faces',meshData.tria(meshData.tnum == 2,1:3),'vertices',meshData.vert, ...
    'facecolor','none', ...
    'edgecolor',[0,0,0]) ;









