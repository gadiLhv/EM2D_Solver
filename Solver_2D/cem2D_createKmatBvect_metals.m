function [K,b] = cem2D_createKmatBvect_metals(meshData,materialList,materialAssign,simProps,f_sim)

% Creates an initial matrix batch using only conducting metals
% Initially only uses surface impedance approximation (see application note ...)
% TBD: Calculate skin depth and use volume resistance instead

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
  
  % If this specific face is dielectric\ferrimagnetic, continue
  if(strcmp(cMaterialProps.type,'metal'))
    continue;
  end
  
  
  % Convert material conductance to project units.
  sigma = cMaterialProps.cond_e;

  % Recover all triangles relevant to current face
  triBinMap = meshData.tnum == faceIdx;
  % Store all triplets 
  triTriplets = meshData.tria(triBinMap,:);
  
  
  % Find edge node pairs in this set of triangles.
  edgePairs = cem2D_findEdgeNodes(triTriplets,meshData.etri);
  
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
  
  % Convert to MKS so Zs_norm will be unitless
  f_mks = units(simProps.freqUnits,'Hz',f_sim);
  
  
  % Equation to follow: d(Ez;Hz)_z/dn = j*k0*((mr/eta);(er*eta))*(Ez;Hz)
  % Alpha coefficients are one
  ax = ay = 1;
  
  
  % According to the polarization, build the coefficients.
  switch simProps.polarizationType
    case 'TE'
      Zs_norm = (1+1i)*sqrt(sigma./(2*pi*f_mks./e0));
      gamma = -1i*k0*er*Zs_norm;
    case 'TM'
      Zs_norm = 1./((1+1i)*sqrt(sigma./(2*pi*f_mks./e0)));
      gamma = -1i*k0*mr*Zs_norm;
    case 'TEM'
      error('''TEM'' currently not supported');
  end
  
  
  % Create interpolants for elements
  ls = cem2D_createEdgeInterpolantCoeffs1st(meshData.vert,edgePairs);
  
  % Initialize triangle element sub-matrices
  Ks = zeros([2 2 size(ls,1)]);
  for i = 1:2
    for j = 1:2
        Ks(i,j,:) = (gamma*ls/6).*(1 + (i == j));
    end
  end

  % Add the elements from the local triangle element sub matrices
  for s = 1:size(Ks,3)
      vIdxs = [edgePairs(s,1) edgePairs(s,2)];
      K(vIdxs,vIdxs) = K(vIdxs,vIdxs) + Ks(:,:,s);
  end
  
  
end