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


 
 
% Start parsing possible responses
narg = numel(varargin);

% Check if each attribute comes with a matching value
if ~~mod(narg,2)
  error('Attributes and values should come in pairs');
end
% Break condition
if narg == 0
  return;
end


% For validation of inputs
validParams = fieldnames(meshPropStruct);

% Update all fields
for argIdx = 1:2:(nargin - 1)
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
%  fprintf(1,['meshPropStruct.' validString ' = varargin{argIdx+1};\n']);
%  setfield(meshPropStruct,validString,varargin{argIdx+1});
  eval(['meshPropStruct.' validString ' = varargin{argIdx+1};']);
  
  if(~ischar(value))
    value = sprintf('%f',value);
  end
  % t wrong possible values
end
  
  
  % TBD: Detect wrong possible values
end