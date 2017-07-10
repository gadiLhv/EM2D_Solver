function meshPropStruct = cem2D_createMeshPropsStruct(varargin)

% Set default options
meshPropStruct = struct(...
    'relWLmeshMax',0.33,...           % 
    'boundingBoxAddSpace',0.25,...    % 
    'useFreespaceWLonly',0,...        % Can be 1 or 0
    'algorithmType','delaunay',...    % This is currently the only one, but 'mesh2D' can cope with others
  );

% Start parsing possible responses
narg = numel(varargin);

% Check if each attribute comes with a matching value
if ~~mod(narg,2)
  error('Attributes and values should come in pairs');
end

% For validation of inputs
validParams = fieldnames(meshPropStruct);

end