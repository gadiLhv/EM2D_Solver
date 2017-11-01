function initMesh = mesh2D_generateInitialMesh(polygonList,materialAssignment,materialProps,meshProps,simProps)

[node,edge,face] = mod2D_polygonToFaceList(polygonList,meshProps.stitchingTolerance);
maxEdgeL = cem2D_calcMaxEdgeL(polygonList{1});

% Define properties for initial lfs (Crude mesh, no sub-meshing)
lfsOpts.kind = meshProps.algorithmType;
lfsOpts.dhdx = meshProps.aspectRatioGrad;
lfsOpts.rho2 = sqrt(maxEdgeL*10);

% Estimate Local Feature Size (LFS) for each node\part. Performing 
% this with 'maxEdgeL' enforces no additional nodes at this point.
[vlfs,tlfs,hlfs] = lfshfn2( node,...
                            edge,...
                            face,...
                            lfsOpts);
                            
% Initial input for the meshing function
lfsStruct = struct('vlfs',vlfs,'tlfs',tlfs,'hlfs',hlfs);

% Initial input for the meshing function
lfsStruct = mesh2D_assignLfsMeshRules(...
              node,...
              edge,...
              face,...
              materialAssignment,...
              lfsStruct,...
              materialProps,...
              meshProps,...
              simProps);


% Super special indexing function for AABB tree queries
slfs = idxtri2(lfsStruct.vlfs,lfsStruct.tlfs);

% Something to do with prior triangular indexing
hfun = @trihfn2;

vlfs = lfsStruct.vlfs;
tlfs = lfsStruct.tlfs;
hlfs = lfsStruct.hlfs;

[ vert,...    % Vertices of mesh cells
  etri,...    % Edges belonging to face (part) boundaries
  tria,...    % Triangle threesomes (attached to VERT)
  tnum] = ... % Part (face) assignments
  refine2(node, ... % Full node list
          edge, ... % Edge connectivity
          face, ... % Edge->Part assignment
          [], ...   % Options field (RHO not used here. Test!)
          hfun,...  % Cool meshing limiting function
          vlfs,...  % Arguments for this function
          tlfs,...
          slfs,...
          hlfs);
          

% Attach all to single structure          
initMesh.vlfs = vlfs;
initMesh.tlfs = tlfs;
initMesh.hlfs = hlfs;
initMesh.slfs = slfs;
initMesh.etri = etri;
initMesh.tria = tria;
initMesh.tnum = tnum;
initMesh.vert = vert;
initMesh.node = node;
initMesh.edge = edge;
initMesh.face = face;

end