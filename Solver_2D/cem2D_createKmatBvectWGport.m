function [K,b] = cem2D_createKmatBvectWGport(meshData,materialList,materialAssign,simProps,f_sim,wgPortStruct)

distTH = 1e-10;

c0 = physical_constant('speed of light in vacuum');
e0 = physical_constant('electric constant');
m0 = physical_constant('mag. constant');

% Convert light speed to correct units
c0 = units('m/sec',[simProps.lengthUnits '*' simProps.freqUnits],c0);
% Calculate wave number
k0 = 2*pi*f_sim/c0;

% 1. Extract from mesh data the relevant edges and nodes
assignedLine = wgPortStruct.assignedLine;
cLineEdges = meshData.lnum == assignedLine;
cLineNodePairs = meshData.ltri(cLineEdges,:);
cLineNodes = cat(3,...
              meshData.vert(cLineNodePairs(:,1),:),...
              meshData.vert(cLineNodePairs(:,2),:));

% Resolve materials onwhich these edges are present


% 1. Detect on which material the WG port "sits"


% 2. This also required to define a mesh "edge"

% 3. Calculate the amplitude values according to the mode #

% 4. Now define the port node values according to the value




end