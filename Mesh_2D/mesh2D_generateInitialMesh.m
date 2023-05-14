function initMesh = mesh2D_generateInitialMesh(polygonList,polylineList,materialAssignment,materialProps,meshProps,simProps)

    [node,edge,face,line,fixes] = mod2D_polygonToFaceList(polygonList,polylineList,meshProps.stitchingTolerance);

    % Maintain only dielectric faces to mesh.
    isFaceDielectric = @(cTxt) strcmp(getfield(cem2D_getMaterialPropsFromName(cTxt,materialProps),'type'),'normal');
    %isFaceDielectric = @(cTxt) strcmp(cTxt,'normal');
    dielectricFaceBin = cellfun(isFaceDielectric,materialAssignment);
    % Store original
    origFaceIdxs = 1:numel(polygonList);
    origMaterialAssignment = materialAssignment;
    origPolygonList = polygonList;
    origFaceList = face;
    % Generate reduced lists
    faceIdxs = origFaceIdxs(dielectricFaceBin);
    materialAssignment = materialAssignment(dielectricFaceBin);
    polygonList = polygonList(dielectricFaceBin);
    face = face(dielectricFaceBin);


    % Function that assigns the maximum edge length with respect to position
    hfun = mesh2D_generateHfun(...
    polygonList,...
    materialAssignment,...
    materialProps,...
    meshProps,...
    simProps);

    % Mesh algorithm options
    refineOpts = struct('RHO2',meshProps.maxRadiusEdgeRatio,...
    'KIND',meshProps.algorithmType,...
    'SIZ1',meshProps.edgeElementTH,...
    'SIZ2',meshProps.triaElementTH);

    [   vert,...    % Vertices of mesh cells
        etri,...    % Edges belonging to face (part) boundaries
        tria,...    % Triangle threesomes (attached to VERT)
        tnum] = ... % Part (face) assignments
        refine2(node, ...         % Full node list
                edge, ...         % Edge connectivity
                face, ...         % Edge->Part assignment
                refineOpts, ...   % Options field (RHO not used here. Test!)
                hfun);            % Cool meshing limiting function

    % Re-assign all the faces\lines\edges to the original list, in case
    % some faces were not meshed during the process (e.g. Metals)

    % Assign edges to polylines
    [   ltri,...    % Lines belonging to ???
        lnum] ...   % Li
        = mesh2D_assignEdgesToPolylines(...
            etri,...
            vert,...
            line,...
            edge,...
            node,...
            meshProps.stitchingTolerance);

    % Face assignment is with respect to the old indexes
    tnum = faceIdxs(tnum);
    tnum = tnum(:);
    face = origFaceList;

    % Attach all to single structure
    initMesh.etri = etri;
    initMesh.tria = tria;
    initMesh.tnum = tnum;
    initMesh.vert = vert;
    initMesh.node = node;
    initMesh.edge = edge;
    initMesh.face = face;
    initMesh.line = line;
    initMesh.ltri = ltri;
    initMesh.lnum = lnum;

end
