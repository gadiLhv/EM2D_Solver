function [K,b] = cem2D_createKmatBvect_2ndOrderRadCond(meshData,meshProps,radEdgeIdxs,materialList,materialAssign,simProps,f_sim)

% Inputs: 
% 1. meshData - Structure with all of the mesh data
% 2. meshProps - General mesh properties
% 3. radEdgeIdxs - Edge numbers in the original faces - meshData.edge(:,edgeIdx) 
% 4. materialList - Cell array with all the material data structures.
% 5. materialAssign - Material assignements per face (object).
% 6. simProps - General simulation properties
% 7. f_sim - Frequency (in simulation units) of current solution

distTH = meshProps.stitchingTolerance;

c0 = physical_constant('speed of light in vacuum');
e0 = physical_constant('electric constant');
m0 = physical_constant('mag. constant');

% Convert light speed to correct units
c0 = units('m/sec',[simProps.lengthUnits '*' simProps.freqUnits],c0);
% Calculate wave number
k0 = 2*pi*f_sim/c0;

edgePairs = meshData.etri(radEdgeIdxs,:);

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
% 4. Extract only outer edge pairscem2D_findEdgeNodes
edgePairs = edgePairs(isOnEdge,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debug: paint triangles and faces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  figure;
%  patch('faces',triTriplets,...
%        'vertices',meshData.vert, ...
%        'facecolor',[1.,1.,1.], ...
%        'edgecolor',[0,0,0]) ;
%  hold on; 
%  axis image off;
%  edgesDir = meshData.vert(edgePairs(:,2),:) - meshData.vert(edgePairs(:,1),:);
%  edgesStart = meshData.vert(edgePairs(:,1),:);  
%  edgesEnd = meshData.vert(edgePairs(:,2),:);  
%  
%  for ei = 1:size(edgesStart,1)
%    quiver(edgesStart(ei,1),edgesStart(ei,2),edgesDir(ei,1)*0.5,edgesDir(ei,2)*0.5,'-b');
%  end
%  hold off;           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debug: paint triangles and faces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate absoulte distance from 
rho1 = meshData.vert(edgePairs(:,1),:);
rho2 = meshData.vert(edgePairs(:,2),:);

% Calculate Kappa coefficients for 'gamma1\2'
kappa = 1./[sqrt(sum(rho1.^2,2)) sqrt(sum(rho2.^2,2))];

% Calculate average gamma values
gamma1 = 1i*k0 + kappa/2 - 1i*(kappa.^2)./(8*(1i*kappa - k0));
gamma2 = -1i./(2*(1i*kappa - k0));
  
gamma1 = 0.5*sum(gamma1,2);
gamma2 = 0.5*sum(gamma2,2);

% Calculate segment lengths
l_s = sqrt(sum((rho2 - rho1).^2,2));

% Construct sub matrice2s
Ks = cat(3,[(gamma1.*l_s/3 + gamma2./l_s) (gamma1.*l_s/6 - gamma2./l_s)],...
           [(gamma1.*l_s/6 - gamma2./l_s) (gamma1.*l_s/3 + gamma2./l_s)]);


% Transform both matrices to matrix batches
Ks = permute(Ks,[3 2 1]);


% Add the elements from the local Ks matrics
for edgeIdx = 1:size(Ks,3)
    vertIdx = edgePairs(edgeIdx,:);
    K(vertIdx,vertIdx) = K(vertIdx,vertIdx) + Ks(:,:,edgeIdx);
end


end

function binAns = inTH(a,b,TH)
binAns = abs(a - b) < TH;
end