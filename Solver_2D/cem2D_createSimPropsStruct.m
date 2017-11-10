function simPropStruct = cem2D_createSimPropsStruct(varargin)

% Set default options
simPropStruct = struct(...
   'fMin',              1 ...               % If nothing else is stated, somethings gotta give
  ,'fMax',              2 ...               % 
  ,'lengthUnits',       'mm' ...            % 
  ,'freqUnits',         'GHz' ...           %  
  ,'timeUnits',         'nsec' ...          %
  ,'numFreqSamps',      1001 ...            %
  ,'adaptiveFreq',      1.5 ...             %
  ,'fieldMonitorFreqs', [] ...              %
  ,'meshProps',         mesh2D_createMeshPropsStruct ...
  ,'polarizationType',  'TEM' ...           % 'TE', 'TM' and 'TEM'.
  ,'solutionType',      'Port Infusion' ... % Can be 'Port Infusion','Eigenmode','Plane Wave'
  ,'minNumOfPasses',    3 ...               % 
  ,'maxNumOfPasses',    6 ...               %
);

% Start parsing possible responses
narg = numel(varargin);

% Check if each attribute comes with a matching value
if ~~mod(narg,2)
  error('Attributes and values should come in pairs');
end

% For validation of inputs
validParams = fieldnames(simPropStruct);

% Update all fields
for argIdx = 1:2:((nargin/2)+1)
  % Validate that paramter has the correct name
  validString = validatestring(varargin{argIdx},validParams);
  % Validate value class
  value = varargin{argIdx+1};
  requiredClass = class(getfield(simPropStruct,validString));
  givenClass = class(value);
  if(~strcmp(requiredClass,givenClass))
    error(...
      sprintf('Parameter ''%s'' needs to be of class ''%s''',validString,requiredClass)...
    );
  end
  
  if(~ischar(value))
    value = sprintf('%f',value);
  end
  
  setfield(simPropStruct,validString,varargin{argIdx+1});
  
  
end

end