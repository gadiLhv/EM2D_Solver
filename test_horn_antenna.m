% test_horn_antenna.m

clc
clear
close all

modelerPath = './Modeler_2D';
mesherPath = '/home/gadi/Repositories/mesh2d';
meshWrapperPath = './Mesh_2D';
solver2DPath = './Solver_2D';
solver1DPath = './Solver_1D';
miscPath = './misc';

addpath(modelerPath);
addpath(genpath(mesherPath));
addpath(meshWrapperPath);
addpath(solver2DPath);
addpath(solver1DPath);
addpath(miscPath);

% Simulation properties
simProps = cem2D_createSimPropsStruct(...
                'fMin',25,...
                'fMax',35,...
                'polarizationType','TE');

f_sim = 30;
                
% Global mesh properties
meshProps = mesh2D_createMeshPropsStruct(...
              'relWLmeshMax',0.33, ...
              'boundingBoxAddSpace',1, ...
              'algorithmType','delfront',...
              'maxRadiusEdgeRatio',1.25);

% Horn antenna properties
horn_WG_L = 10;
horn_WG_H = 5;
horn_taper_ang_deg = 30;
horn_taper_L = 6;
horn_metal_thick = 0.8;

%%%%%%%%%%%%%%%%%%%%%%%
% Build horn geometry %
%%%%%%%%%%%%%%%%%%%%%%%

% Draw top and bottom panels
topPanel = mod2D_createRectangleStruct([-0.5*horn_WG_L 0.5*horn_WG_H],[0.5*horn_WG_L 0.5*horn_WG_H+horn_metal_thick]);
botPanel = mod2D_createRectangleStruct([-0.5*horn_WG_L -0.5*horn_WG_H-horn_metal_thick],[0.5*horn_WG_L -0.5*horn_WG_H]);

% 
backPlate = mod2D_createRectangleStruct([0.5*horn_WG_L -0.5*horn_WG_H-horn_metal_thick],[0.5*horn_WG_L+horn_metal_thick 0.5*horn_WG_H+horn_metal_thick]);
portLine = mod2D_createLineStruct([0.5*horn_WG_L -0.5*horn_WG_H-horn_metal_thick],[0.5*horn_WG_L 0.5*horn_WG_H+horn_metal_thick]);

taper_ang = horn_taper_ang_deg*pi/180;

topTaper = mod2D_createPolygonStruct(...
            [-0.5*horn_WG_L ; ...
             -0.5*horn_WG_L-horn_taper_L ; ...
             -0.5*horn_WG_L-horn_taper_L ; ...
             -0.5*horn_WG_L], ...
            [0.5*horn_WG_H ; ...
             0.5*horn_WG_H+horn_taper_L*tan(taper_ang) ; ...
             0.5*horn_WG_H+horn_taper_L*tan(taper_ang)+horn_metal_thick ; ...
             0.5*horn_WG_H+horn_metal_thick]);

botTaper = mod2D_createPolygonStruct(...
            [-0.5*horn_WG_L ; ...
             -0.5*horn_WG_L-horn_taper_L ; ...
             -0.5*horn_WG_L-horn_taper_L ; ...
             -0.5*horn_WG_L], ...
            [-0.5*horn_WG_H ; ...
             -0.5*horn_WG_H-horn_taper_L*tan(taper_ang) ; ...
             -0.5*horn_WG_H-horn_taper_L*tan(taper_ang)-horn_metal_thick ; ...
             -0.5*horn_WG_H-horn_metal_thick]);
% Show geometry
topPanel = mod2D_booleanOperation(topPanel,topTaper,'add');
botPanel = mod2D_booleanOperation(botPanel,botTaper,'add');

%%%%%%%%%%%%%%%%%%%%
% Assign materials %
%%%%%%%%%%%%%%%%%%%%

% Create the materials for this simulation
material_default = cem2D_createMaterialDefs;
material_PEC = cem2D_createMaterialDefs('name','PEC','type','PEC');
% Create the material list with only the default material
materialList = cem2D_addMaterialToList(material_default);
materialList = cem2D_addMaterialToList(material_PEC,materialList);

% Create polygon list 
polList = {topPanel,botPanel,backPlate};
% Assign materials per-shape
materialAssignement = {'PEC','PEC','PEC'};

% Create bounding box (background)
bBox = mod2D_createBoundingBox(...
        polList,...
        meshProps,...
        simProps,...
        material_default);
        
% Subtract all shapes from bounding box
for polIdx = 1:numel(polList)
    bBox = mod2D_booleanOperation(bBox,polList{polIdx},'subtract');
end

polList = [{bBox} polList];
lineList = {portLine};
materialAssignement = [{'default'},materialAssignement];


%%%%%%%%%%%%%%%%%%%%%
% Plot the geometry %
%%%%%%%%%%%%%%%%%%%%%
figHdl = figure;
axHdl = axes;
mod2D_showPolygon(axHdl,bBox,[1 1 1],[0 0 0]);
mod2D_showPolygon(axHdl,topPanel,[0.3 0.3 0.3],[0 0 0]);
mod2D_showPolygon(axHdl,botPanel,[0.3 0.3 0.3],[0 0 0]);
mod2D_showPolygon(axHdl,backPlate,[0 0.1 0.7],[0 0 0]);

set(gca,'fontsize',14);
xlabel('x [mm]','fontsize',16);
ylabel('y [mm]','fontsize',16);
%%%%%%%%%%%%%%%%%%%%%
% Plot the geometry %
%%%%%%%%%%%%%%%%%%%%%

% Test mesh of bounding box
meshData = mesh2D_generateInitialMesh(...
            polList,...             % Entire polygon list. This is mainly to calculate the maximum\minimum edge length
            lineList,...            % List of lines. Used for porst, mainly.
            materialAssignement,... % Material assignments for initial LFS assignment
            materialList,...        % Corresponding material properties
            meshProps,...           % Darrens Mesher properties
            simProps);              % Simulation properties. This is to assign spatial meshing rules, most

% Smooth the mesh
meshData = mesh2D_smoothMesh(meshData,meshProps);

hold(axHdl,'on');
patch('faces',meshData.tria(meshData.tnum == 1,1:3),'vertices',meshData.vert, ...
    'facecolor','none', ...
    'edgecolor',[0,0,0]) ;

% Plot lines
lineIdxs = unique(meshData.lnum);
for cLineIdx = lineIdxs(:).'
    cLineVertPairs = meshData.ltri(meshData.lnum == cLineIdx,:);
    
    cLineR1 = meshData.vert(cLineVertPairs(:,1),:);
    cLineR2 = meshData.vert(cLineVertPairs(:,2),:);
    
    plot([cLineR1(:,1) cLineR2(:,1)].',[cLineR1(:,2) cLineR2(:,2)].','-r','linewidth',3);
    plot(cLineR1(:,1),cLineR1(:,2),'.b','markersize',20);
    plot(cLineR2(:,1),cLineR2(:,2),'.g','markersize',10);
end

hold(axHdl,'off');

%mod2D_showPolyline(axHdl,portLine,[1 0 0]);    

portStruct = cem2D_createWaveguidePortStruct('polyLineNumber',1);

% For this example  , set the simulation frequency as the middle between
% minimum and maximum frequency
simProps.fSim = 0.5*(simProps.fMin + simProps.fMax);

% Generate port modes, store in port structure.
portStruct = cem1D_calcPortModes(...
                portStruct,...
                meshData,...
                materialList,...
                materialAssignement,...
                meshProps,...
                simProps);

% Switch propagation direction
portStruct.propDirection = -1;

% Show propagation direction
hold(axHdl,'on');
arrowHdl = quiver(  portStruct.portCenter(1),...
                    portStruct.portCenter(2),...
                    4*portStruct.portDirection(1)*portStruct.propDirection,...
                    4*portStruct.portDirection(2)*portStruct.propDirection);
                    

set(arrowHdl,'linewidth',3);
set(arrowHdl,'color',[1 0.2 0.2]);
hold(axHdl,'off');

% Create matrix for port propagation
[K_port,b_port] = cem2D_createKmatBvect_WG_Port(...
                    meshData,...
                    portStruct,...
                    materialList,...
                    simProps,...
                    meshProps,...
                    f_sim);


%% Assign to outline the 2nd order radiating boundary conditions
%[K_rad,b_rad] = cem2D_createKmatBvect_2ndOrderRadCond(...
%                    smoothMesh,...              % Mesh to use for final matrix
%                    materialList,...            % List of all materials used
%                    materialAssignement,...     % Material assignements per-face
%                    simProps,...                % Simulation properties
%                    meshProps,...               % Mesh properties
%                    0.5*(simProps.fMax + simProps.fMin));
%                    
rmpath(modelerPath);
rmpath(genpath(mesherPath));
rmpath(meshWrapperPath);
rmpath(solver2DPath);
rmpath(solver1DPath);
rmpath(miscPath);