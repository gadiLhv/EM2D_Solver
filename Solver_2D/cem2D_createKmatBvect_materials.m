function [K,b] = cem2D_createKmatBvect_materials(meshData,materialList,materialAssign,simProps,f_sim)
% Creates an initial matrix batch using only dielectric and ferrimagnetic
% properties of the medium. This formulation is relevant for TM(Z) polarization
%
% Inputs: 
% 1. meshData - Structure with all of the mesh data
% 2. materialList - Cell array with all the material data structures.
% 3. materialAssign - Material assignements per face (object).
% 4. simProps - General simulation properties
% 5. f_sim - Frequency (in simulation units) of current solution

c0 = physical_constant('speed of light in vacuum');
e0 = physical_constant('electric constant');
m0 = physical_constant('mag. constant');

% Convert light speed to correct units
c0 = units('m/sec',[simProps.lengthUnits '*' simProps.freqUnits],c0);
% Calculate wave number
k0 = 2*pi*f_sim/c0;

nVerts = size(meshData.vert,1);

% Initial FEM matrix and source vector
K = zeros([1 1]*nVerts);
b = zeros([nVerts 1]);
  
% Build matrix face by face
for faceIdx = 1:numel(meshData.face)
  % Extract current parts material properties
  cMaterialProps = cem2D_getMaterialPropsFromName(materialAssign{faceIdx},materialList);
  
  
  % If this specific face is metal, continue
  if(strcmp(cMaterialProps.type,'metal'))
    continue;
  end
  
  % Get relative permiability and permittivity
  mr = cMaterialProps.mr;
  er = cMaterialProps.er;
  
  % Recover all triangles relevant to current face
  triBinMap = meshData.tnum == faceIdx;
  % Store all triplets 
  triTriplets = meshData.tria(triBinMap,:);
  
  % Find edge node pairs in this set of 
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
%  hold off;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Debug: paint triangles and faces %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%triangles.
  edgePairs = cem2D_findEdgeNodes(triTriplets,meshData.etri);
  cem2D_createKmatBvect_2ndOrderRadCond
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
%  hold off;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Debug: paint triangles and faces %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % According to the polarization, build the coefficients.
  switch simProps.polarizationType
    case 'TE'
      ax = ay = 1/er;
      beta = mr*k0^2;
    case 'TM'
      ax = ay = 1/mr;
      beta = er*k0^2;
    case 'TEM'
      error('''TEM'' currently not supported');
  end
  
  % Create interpolants for elements
  [a,b,c,Det] = cem2D_createInterpolantCoeffs1st(meshData.vert,triTriplets);
  
  % Initialize triangle element sub-matrices
  Ke = zeros([3 3 size(c,1)]);
  for i = 1:3
    for j = 1:3
        Ke(i,j,:) = ((4*Det).^(-1)).*(ax.*b(:,i).*b(:,j) + ...
                    ay.*c(:,i).*c(:,j)) + (Det/12).*beta.*(1 + (i == j));
    end
  end

  % Add the elements from the local triangle element sub matrices
  for e = 1:size(ke,3)
      vIdxs = [triTriplets(e,1) triTriplets(e,2) triTriplets(e,3)];
      K(vIdxs,vIdxs) = K(vIdxs,vIdxs) + Ke(:,:,e);
  end
  
  
end



end