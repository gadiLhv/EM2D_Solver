function lfsStruct = cem2D_assignLfsMeshRules(polygonList,materialAssignment,lfsStruct,materialList,meshProps,simProps)
% After Local Feature Size (LFS) approximation, there is a new triangular 
% assignement. According to material properties, lfsh rules are to be changed,
% according to "worse case" (shortest wavelength).

c0 = physical_constant('speed of light in vacuum');
e0 = physical_constant('electric constant');
m0 = physical_constant('mag. constant');

vlfs = lfsStruct.vlfs;    % Vertices of initial LFS
tlfs = lfsStruct.tlfs;    % Triangles constructed by initial LFS
hlfs = lfsStruct.hlfs;    % LFS approximated by initial LFS


% Preliminary data to determine mesh sizes
c0 = meshProps.c0;    % Free space speed of light
f0 = meshProps.f0;    % Frequency for current mesh

% Extract threshold for "Same node"\"Node on edge" criteria
distTH = meshProps.sameNodeDistTH;

% Relative mesh edge maximum with respect to wavelength
relWLmeshMax = meshProps.relWLmeshMax; 

wl0 = c0/f0;          % Free space wavelength.

% Normalize to units
wl0 = wl0*mod2D_lengthUnitFactor(meshProps.lengthUnits);

% Extract relative permittivity and permiability and assign the maximum 
% allowed mesh size.
edgeMax = [];
for faceIdx = 1:numel(face)
  % Extract current parts material properties
  cMaterialProps = CEM2D_getMaterialPropsFromName(materialAss{faceIdx},materialList);
  m0 = cMaterialProps.m0;
  e0 = cMaterialProps.e0;
  
  % Set maximum wavelength fraction of this face's mesh.
  edgeMax = [edgeMax ; wl0*relWLmeshMax/sqrt(e0*m0)];
end

% For each face, find all the relevant nodes
for faceIdx = 1:numel(face)
  % Start from geometry node list
  cFace = face{faceIdx};
  cEdge = edge(cFace,:);
  
  % Now match the vertices from the lfs list to this face
  [vertIdx,~,~] = mod2D_pointOnEdge(vlfs,...      % All the vertices in the LFS generation
                                    node,...      % Full node list in the original geometry
                                    cEdge,...     
                                    distTH);
  
  % For each of these vertices, match a new maximum edge length.
  cHlfs = hlfs(vertIdx);
  cHlfs = min(cHlfs,edgeMax(faceIdx));
  hlfs(vertIdx) = cHlfs;
  
end

% Update the lfs structure
lfsStruct.hlfs = hlfs;

end