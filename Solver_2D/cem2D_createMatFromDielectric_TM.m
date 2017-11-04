function [K,b] = cem2D_createKmat_TM(meshData,materialList,materialAssign,simProps,f0)
% Creates an initial matrix batch using only dielectric and ferrimagnetic
% properties of the medium. This formulation is relevant for TM(Z) polarization
%
% Inputs: 
% 1. meshData - Structure with all of the mesh data
% 2. materialList - Cell array with all the material data structures.
% 3. materialAssign - Material assignements per face (object).
% 4. simProps - General simulation properties

c0 = physical_constant('speed of light in vacuum');
e0 = physical_constant('electric constant');
m0 = physical_constant('mag. constant');

edgeCouples = meshData.etri;

% Build matrix face by face
for faceIdx = 1:numel(meshData.face)
  % Extract current parts material properties
  cMaterialProps = cem2D_getMaterialPropsFromName(materialAssign{faceIdx},materialList);
  
  % Get relative permiability and permittivity
  mr = cMaterialProps.mr;
  er = cMaterialProps.er;
  
  % Recover all triangles relevant to current face
  triIdxs = meshData.tnum == faceIdx;
  % Store all triplets 
  triTriplets = meshData.tria(triIdxs,:);
  
  % Recover all edge nodes corresponding to the
  % extracted triangles:
  % 1. Determine if either of the edge nodes is present in a triangle
  edgeNode1map = bsxfun(@eq,edgeCouples(:,1),permute(triTriplets,[3 1 2]));
  edgeNode2map = bsxfun(@eq,edgeCouples(:,2),permute(triTriplets,[3 1 2]));
  % 2. Determine if both of them are present in a triangle
  areBothNodesIn = (sum(edgeNode1map,3) + sum(edgeNode2map,3)) == 2;
  % 3. Extract only nodes that consist an edge in this face.
  
  
end


end