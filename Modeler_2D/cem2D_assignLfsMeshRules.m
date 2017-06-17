function lfsStruct = cem2D_assignLfsMeshRules(geomStruct,lfsStruct,materialProps,meshProps)
% After Local Feature Size (LFS) approximation, there is a new triangular 
% assignement. According to material properties, lfsh rules are to be changed,
% according to "worse case" (shortest wavelength).

node = geomStruct.node;   % Current geometry initial node list
edge = geomStruct.edge;   % Current geometry initial edge list
face = geomStruct.face;   % Current geometry edge->face assignement cell array

vlfs = lfsStruct.vlfs;    % Vertices of initial LFS
tlfs = lfsStruct.tlfs;    % Triangles constructed by initial LFS
hlfs = lfsStruct.hlfs;    % LFS approximated by initial LFS

% Preliminary data to determine mesh sizes
c0 = meshProps.c0;    % Free space speed of light
f0 = meshProps.f0;    % Frequency for current mesh

% Relative mesh edge maximum with respect to wavelength
relWLmeshMax = meshProps.relWLmeshMax; 

wl0 = c0/f0;          % Free space wavelength.

% Normalize to units
wl0 = wl0*mod2D_lengthUnitFactor(meshProps.lengthUnits);

% Extract relative permittivity and permiability and assign the maximum 
% allowed mesh size.
edgeMax = [];
for faceIdx = 1:numel(face)
  cMaterialProps = materialProps{faceIdx};
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
  [vertIdx,whichEdge,t] = mod2D_pointOnEdge(vlfs,node,edge,1e-9);
  
  % For each of these vertices, match a new maximum edge length.
  cHlfs = hlfs(vertIdx);
  cHlfs = min(cHlfs,edgeMax(faceIdx));
  hlfs(vertIdx) = cHlfs;
  
end

% Update the lfs structure
lfsStruct.hlfs = hlfs;

end