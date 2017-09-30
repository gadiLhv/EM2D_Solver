function meshPropStruct = cem2D_createMeshPropsStruct(varargin)

% Set default options
meshPropStruct = struct(...
    'relWLmeshMax',0.33,...           % 
    'boundingBoxAddSpace',0.25,...    % 
    'useFreespaceWLonly',0,...        % Can be 1 or 0
    'algorithmType','delaunay',...    % This is currently the only one, but 'mesh2D' can cope with others
    'stitchingTolerance',1e-8, ...    % 
    'performMeshSmoothing',1, ...     % 
  );


 
 
% Start parsing possible responses
narg = numel(varargin);

% Check if each attribute comes with a matching value
if ~~mod(narg,2)
  error('Attributes and values should come in pairs');
end

% For validation of inputs
validParams = fieldnames(meshPropStruct);

% Update all fields
for argIdx = 1:2:((nargin/2)+1)
  % Validate that paramter has the correct name
  validString = validatestring(varargin{argIdx},validParams);
  % Validate value class
  value = varargin{argIdx+1};
  requiredClass = class(getfield(meshPropStruct,validString));
  givenClass = class(value);
  if(~strcmp(requiredClass,givenClass))
    error(...
      sprintf('Parameter ''%s'' needs to be of class ''%s''',validString,requiredClass)...
    );
  end
  
  setfield(meshPropStruct,validString,varargin{argIdx+1});
  
  if(~ischar(value))
    value = sprintf('%f',value);
  end
  t wrong possible values
end
  eval(['meshPropStruct. ' validString ' = ' value ';']);
  
  % TBD: Detect wrong possible values
end

end