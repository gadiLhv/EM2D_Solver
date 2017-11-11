function [K,b] = cem2D_createKmatBvect_2ndOrderRadCond(meshData,materialList,materialAssign,simProps,f_sim)

% Inputs: 
% 1. meshData - Structure with all of the mesh data
% 2. materialList - Cell array with all the material data structures.
% 3. materialAssign - Material assignements per face (object).
% 4. simProps - General simulation properties
% 5. f_sim - Frequency (in simulation units) of current solution

distTH = 1e-10;

c0 = physical_constant('speed of light in vacuum');
e0 = physical_constant('electric constant');
m0 = physical_constant('mag. constant');

% Convert light speed to correct units
c0 = units('m/sec',[simProps.lengthUnits '*' simProps.freqUnits],c0);
% Calculate wave number
k0 = 2*pi*f_sim/c0;

% Extract background material (allways first face) properties
cMaterialProps = cem2D_getMaterialPropsFromName(materialAssign{1},materialList);

nVerts = size(meshData.vert,1);

% Initial FEM matrix and source vector
K = zeros([1 1]*nVerts);
b = zeros([nVerts 1]);
  

% Get relative permiability and permittivity
mr = cMaterialProps.mr;
er = cMaterialProps.er;

% Recover all triangles relevant to background face
triBinMap = meshData.tnum == 1;
% Store all triplets 
triTriplets = meshData.tria(triBinMap,:);

% Find edge node pairs in this set of triangles.
edgePairs = cem2D_findEdgeNodes(triTriplets,meshData.etri);

% Extract only edges on the outermost part
% 1. Find the extrema of the bounding box
xRange = [min(meshData.vert(:,1)) max(meshData.vert(:,1))];
yRange = [min(meshData.vert(:,2)) max(meshData.vert(:,2))];
% 2. Extract all pairs
Rpairs1 = meshData.vert(edgePairs(:,1),:);
Rpairs2 = meshData.vert(edgePairs(:,2),:);
% 3. Determine if the pairs are on any of the extemum. Decide with threshold, 
%    rather than simply comparing. Mesh isn't perfect...
isOnEdge =  (inTH(Rpairs1(:,1),xRange(1),distTH) & inTH(Rpairs2(:,1),xRange(1),distTH)) | ...
            (inTH(Rpairs1(:,1),xRange(2),distTH) & inTH(Rpairs2(:,1),xRange(2),distTH)) | ...
            (inTH(Rpairs1(:,2),yRange(1),distTH) & inTH(Rpairs2(:,2),yRange(1),distTH)) | ...
            (inTH(Rpairs1(:,2),yRange(2),distTH) & inTH(Rpairs2(:,2),yRange(2),distTH));
% 4. Extract only outer edge pairs
edgePairs = edgePairs(isOnEdge,:);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Debug: paint triangles and faces %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  figure;
%  patch('faces',triTriplets,...
%        'vertices',meshData.vert, ...
%        'facecolor',[1.,1.,1.], ...
%        'edgecolor',[0,0,0]) ;
%  hold on; 
%  axis image off;
%  edgesX = [meshData.vert(edgePairs(:,1),1) meshData.vert(edgePairs(:,2),1)].';
%  edgesY = [meshData.vert(edgePairs(:,1),2) meshData.vert(edgePairs(:,2),2)].';
%  plot(edgesX,edgesY,'-','linewidth',3);
%  hold off;       'facecolor',[1.,1.,1.], ...
        'edgecolor',[0,0,0]) ;
  hold on; 
  axis image off;
  edgesX = [meshData.vert(edgePairs(:,1),1) meshData.vert(edgePairs(:,2),1)].';
  edgesY = [meshData.vert(edgePairs(:,1),2) meshData.vert(edgePairs(:,2),2)].';
  plot(edgesX,edgesY,'-','linewidth',3);
  hold off;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Debug: paint triangles and faces %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%% Create interpolants for elements
%[a,b,c,Det] = cem2D_createInterpolantCoeffs1st(meshData.vert,triTriplets);



end

function binAns = inTH(a,b,TH)
binAns = abs(a - b) < TH;
end