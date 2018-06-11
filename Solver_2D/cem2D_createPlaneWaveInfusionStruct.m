function pwStruct = cem2D_createPlaneWaveInfusionStruct(varargin)
% Polarization is set via simulation properties

pwStruct = struct('phiPropagation',0,...  % Phi of propagation direction, in degrees
                  'amp',1,...             % Amplitude in V/m (for E-field. Need to be corrected for 0 field
                  'polarizationVect',[1 0],...  % Computed around propagatin axis, in transverse to the computation plane
                  'measurementPlaneDist',0,...  % Phase is zero at this plane. Defaulted in [0,0].
                  'measurementStartingPoint',[0 0],... % Distance is measured from this coordinate
                  );
    

% For validation of inputs
validParams = fieldnames(pwStruct);

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
  if(~ischar(value))
    value = sprintf('%f',value);
  end
  
  setfield(simPropStruct,validString,varargin{argIdx+1});
  
  
end
    
end