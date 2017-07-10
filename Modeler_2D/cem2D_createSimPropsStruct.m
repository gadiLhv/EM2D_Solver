function simPropStruct = cem2D_createSimPropsStruct(varargin)

% Set default options
simPropStruct = struct(...
  'fMin',1,...              % If nothing else is stated, somethings gotta give
  'fMin',2,...              % 
  'lengthUnits','mm',...    % 
  'freqUnits','GHz',...     
  'timeUnits','nsec',...
);

% Start parsing possible responses
narg = numel(varargin);

% Check if each attribute comes with a matching value
if ~~mod(narg,2)
  error('Attributes and values should come in pairs');
end

% For validation of inputs
validParams = fieldnames(simPropStruct);



end