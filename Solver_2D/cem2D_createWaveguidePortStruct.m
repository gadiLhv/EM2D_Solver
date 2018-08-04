function wgPortStruct = cem2D_createWaveguidePortStruct(varargin)
  
wgPortStruct = struct('nMode',0,...                 % Varies number of peaks. 0 for discrete port, 1 for single peak sinus
                      'power',1,...                 % Port power. Normalized to 50Ohm
                      'measurementPlaneDist',0,...  % De-embedding
                      'normImpedance',50,...        % Impedance to re-normalize
                      'assignedLine',1              % Index of the assigned polyline line assignement list.
                  );
                  

% For validation of inputs
validParams = fieldnames(wgPortStruct);

% Update all fields
for argIdx = 1:2:((nargin/2)+1)
  % Validate that paramter has the correct name
  validString = validatestring(varargin{argIdx},validParams);
  % Validate value class
  value = varargin{argIdx+1};
  requiredClass = class(getfield(pwStruct,validString));
  givenClass = class(value);
  if(~strcmp(requiredClass,givenClass))
    error(...
      sprintf('Parameter ''%s'' needs to be of class ''%s''',validString,requiredClass)...
    );
  end
  
  if(~ischar(value))
    value = sprintf('%f',value);
  end
  
  setfield(wgPortStruct,validString,varargin{argIdx+1});

end