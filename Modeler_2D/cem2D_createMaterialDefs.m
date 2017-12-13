function materialDef = cem2D_createMaterialDefs(varargin)

% Set default options
materialDef = struct(...
  'name','default',...  % Material name
  'er',1,...            % Relative permittivity
  'mr',1,...            % Relative permiability
  'type','normal',...   % Could be 'normal' or 'metal'. Namely: Dielectric or metallic
  'cond_e',1e20,...     % Electric conductivity - used only in case of metallic. Given in Ohm\m
  'cond_m',1e20,...     % Magnetic conductivity - used only in case of metallic. Given in ??? (TBD)
  'tand_e',0,...        % Dielectric loss tangent - used only in case of normal
  'tand_m',0,...        % Magnetic loss tangent - used only in case of normal
  'f_m',0 ...           % The frequency in which the data is measured
  );

% Start parsing possible responses
narg = numel(varargin);

% Check if each attribute comes with a matching value
if ~~mod(narg,2)
  error('Attributes and values should come in pairs');
end

% For validation of inputs
validParams = fieldnames(materialDef);

% Iterate and replace values
for argIdx = 1:(numel(varargin)/2)
  % Extract name and value
  argName = varargin{argIdx*2-1};
  argVal = varargin{argIdx*2};
  
  % Validate the attribute
  argName = validatestring(argName,validParams,'cem2D_createMaterialDefs','attribute');
  
  % Parse inputs
  switch argName
    case 'name'
      if ~ischar(argVal)
        error(sprintf('Matching value to argument ''%s'' needs to be a character string',argName));
      end
      materialDef.name = argVal;
    case 'type'
      % Valudate 'normal' or 'metal'
      argVal = validatestring(argVal,{'normal','metal'},'cem2D_createMaterialDefs','type');
      materialDef.type = tolower(argVal);
    otherwise
      % Validate numeric
      if ~isnumeric(argVal)
        error(sprintf('Matching value to argument ''%s'' needs to be numeric',argName));
      end
      % Assign value
      strToEval = sprintf('materialDef.%s = %f;',argName,argVal);
      eval(strToEval);
  end
end


end
