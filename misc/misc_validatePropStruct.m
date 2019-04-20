function propStruct = misc_validatePropStruct(propStruct,arglist)
    narg = numel(arglist);
    
        % Check if each attribute comes with a matching value
    if ~~mod(narg,2)
      error('Attributes and values should come in pairs');
    end

    % For validation of inputs
    validParams = fieldnames(propStruct);

    % Update all fields
    for argIdx = 1:2:(narg-1)
      % Validate that paramter has the correct name
      validString = validatestring(arglist{argIdx},validParams);
      % Validate value class
      value = arglist{argIdx+1};
      requiredClass = class(getfield(propStruct,validString));
      givenClass = class(value);
      if(~strcmp(requiredClass,givenClass))
        error(...
          sprintf('Parameter ''%s'' needs to be of class ''%s''',validString,requiredClass)...
        );
      end
      
      eval(['propStruct.' validString ' = arglist{argIdx+1};']);
      
    end

end
