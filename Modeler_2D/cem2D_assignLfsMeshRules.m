function lfsStruct = cem2D_assignLfsMeshRules(node,edge,face,materialAssignment,lfsStruct,materialList,meshProps,simProps)
% After Local Feature Size (LFS) approximation, there is a new triangular 
% assignement. According to material properties, lfsh rules are to be changed,
% according to "worse case" (shortest wavelength).

c0 = physical_constant('speed of light in vacuum');
e0 = physical_constant('electric constant');
m0 = physical_constant('mag. constant');

% TBD: Currently mesh size is set according to the highest simulated frequency.
f0 = simProps.fMax;

vlfs = lfsStruct.vlfs;    % Vertices of initial LFS
tlfs = lfsStruct.tlfs;    % Triangles constructed by initial LFS
hlfs = lfsStruct.hlfs;    % LFS approximated by initial LFS

% Extract threshold for "Same node"\"Node on edge" criteria
distTH = meshProps.stitchingTolerance;

% Relative mesh edge maximum with respect to wavelength
relWLmeshMax = meshProps.relWLmeshMax; 

% Free space wavelength.
c0 = units('m/sec',[simProps.lengthUnits '*' simProps.freqUnits],c0);

% Wavelength with respect to current units
wl0 = c0/f0;

% Extract relative permittivity and permiability and assign the maximum 
% allowed mesh size.
edgeMax = [];
for faceIdx = 1:numel(face)
  % Extract current parts material properties
  cMaterialProps = cem2D_getMaterialPropsFromName(materialAssignment{faceIdx},materialList);
  
  % Get relative permiability and permittivity
  mr = cMaterialProps.mr;
  er = cMaterialProps.er;
  
  % Set maximum wavelength fraction of this face's mesh.
  edgeMax = [edgeMax ; wl0*relWLmeshMax/sqrt(er*mr)];
end

% For each face, find all the relevant nodes
for faceIdx = 1:numel(face)
  % Start from geometry node list
  cFace = face{faceIdx};
  cEdge = edge(cFace,:);
  
  % Now match the vertices from the lfs list to this face
  [vertIdx,~,~] = mod2D_pointOnEdge(vlfs,...   % All the vertices in the LFS generation
                                    node,...   % Full node list in the original geometry
                                    cEdge,...     
                                    distTH);
  
  % Dont repeat nodes
  vertIdx = unique(vertIdx);
  
  % For each of these vertices, match a new maximum edge length.
  cHlfs = hlfs(vertIdx);
  cHlfs = min(cHlfs,edgeMax(faceIdx));
  hlfs(vertIdx) = cHlfs;
  
end

% Update the lfs structure
lfsStruct.hlfs = hlfs;

end