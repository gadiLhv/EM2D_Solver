function cem2D_createSolverSettingsStruct(varargin)

solverSettingsStruct = struct(...
  'freqRange',[1 2] ...
  ,'numFreqSamps',1001 ...
  ,'adaptiveFreq',1.5 ...
  ,'fieldMonitorFreqs',[] ...
  ,'meshProps',cem2D_createMeshPropsStruct,...
  ,'polarizationType','TEM' ...
  ,'solutionType','Port Infusion' ... % Can be 'Port Infusion' and 'Eigenmode'
  ,'minNumOfPasses',3 ...
  ,'maxNumOfPasses',6 ...
  );
  
% Start parsing possible responses
narg = numel(varargin);

% Check if each attribute comes with a matching value
if ~~mod(narg,2)
  error('Attributes and values should come in pairs');
end

% For validation of inputs
validParams = fieldnames(solverSettingsStruct);

% Update all fields
for argIdx = 1:2:((nargin/2)+1)
  % Validate that paramter has the correct name
  validString = validatestring(varargin{argIdx},validParams);
  % Validate value class
  value = varargin{argIdx+1};
  requiredClass = class(getfield(solverSettingsStruct,validString));
  givenClass = class(value);
  if(~strcmp(requiredClass,givenClass))
    error(...
      sprintf('Parameter ''%s'' needs to be of class ''%s''',validString,requiredClass)...
    );
  end
  
  setfield(solverSettingsStruct,validString,varargin{argIdx+1});
  
end

end
