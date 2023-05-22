% test_microstrip_line.m


% test_horn_antenna.m

warning off;
clc
clear
close all

pkg load geometry
pkg load miscellaneous

modelerPath = './Modeler_2D';
mesherPath = '../mesh2d';
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
                'fMin',0.1,...
                'fMax',5,...
                'polarizationType','TE');

                % Global mesh properties
meshProps = mesh2D_createMeshPropsStruct(...
              'relWLmeshMax',0.1, ...
              'boundingBoxAddSpace',1, ...
              'algorithmType','delfront',...
              'maxRadiusEdgeRatio',1.5);

% Unless overridden, f_sim will be set to the average between fMin and fMax
% f_sim = 3;

% Line Properties

WG_H = 55;
WG_W = 110;

WG_Box = mod2D_createRectangleStruct([-0.5*WG_W 0],[0.5*WG_W WG_H]);

% Create the materials for this simulation
material_default = cem2D_createMaterialDefs;
material_PEC = cem2D_createMaterialDefs('name','PEC','type','PEC');
material_FR4 = cem2D_createMaterialDefs('name','FR4','type','normal','er',4.3);
% Create the material list with only the default material
materialList = cem2D_addMaterialToList(material_default);
materialList = cem2D_addMaterialToList(material_PEC,materialList);
materialList = cem2D_addMaterialToList(material_FR4,materialList);

% Create polygon list
polList = {WG_Box};

% Assign materials per-shape
materialAssignement = {'default'};


%%%%%%%%%%%%%%%%%%%%%
% Plot the geometry %
%%%%%%%%%%%%%%%%%%%%%
figHdl = figure;
axHdl = axes;
mod2D_showPolygon(axHdl,WG_Box,[1 1 1],[0 0 0]);

set(gca,'fontsize',14);
xlabel('x [mm]','fontsize',16);
ylabel('y [mm]','fontsize',16);
%%%%%%%%%%%%%%%%%%%%%
% Plot the geometry %
%%%%%%%%%%%%%%%%%%%%%

lineList = {};

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

if exist('f_sim','var')
    if isempty(f_sim)
        f_sim = 0.5*(simProps.fMin + simProps.fMax);
    end
else
    f_sim = 0.5*(simProps.fMin + simProps.fMax);
end

[kz,fc,Et,Ez,edgeTriplets] = cem2D_calcPortModes(meshData,meshProps,materialList,materialAssignement,simProps,f_sim);

hold(axHdl,'on');

for tIdx = unique(meshData.tnum).'
    patch('faces',meshData.tria(meshData.tnum == tIdx,1:3),'vertices',meshData.vert, ...
        'facecolor','none', ...
        'edgecolor',[0,0,0]) ;
end

[Xm,Ym] = meshgrid(linspace(-0.5*WG_W,0.5*WG_W,25),linspace(0,WG_H,15));

EI = cem2D_vectorElementInterp(...
    meshData.vert,...
    meshData.tria,...
    edgeTriplets,...
    Et(:,1));

Exy = EI(Xm,Ym);

qHdl = quiver(Xm(:),Ym(:),abs(Exy(:,1)),abs(Exy(:,2)));
set(qHdl,'color',[1 0 0]);

hold(axHdl,'off');

warning off;
