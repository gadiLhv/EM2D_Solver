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

    simPropStruct = misc_validatePropStruct(simPropStruct,varargin);

end