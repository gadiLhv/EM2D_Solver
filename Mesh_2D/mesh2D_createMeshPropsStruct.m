function meshPropStruct = mesh2D_createMeshPropsStruct(varargin)

    % Set default options
    meshPropStruct = struct(...
        'relWLmeshMax',0.33,...                 % 
        'boundingBoxAddSpace',0.25,...          % 
        'useFreespaceWLonly',0, ...             % Can be 1 or 0
        'stitchingTolerance',1e-8, ...          % Relevant for geometry pre-processing
        'algorithmType','delfront',...          % 'delauny' or 'delfront' - Relevant for both lfs and refine2 (and smoothing?)
        'maxRadiusEdgeRatio',1.025, ...  % Read in 'refine2' (mesh2D) about RHO2
        'edgeElementTH',1.333, ...       % Read in 'refine2' (mesh2D) about SIZ1
        'triaElementTH',1.3, ...         % Read in 'refine2' (mesh2D) about SIZ1
        'performMeshSmoothing',1  ...           % 
      );


 
    meshPropStruct = misc_validatePropStruct(meshPropStruct,varargin);

end