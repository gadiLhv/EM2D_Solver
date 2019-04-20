function portStruct = cem2D_createWaveguidePortStruct(varargin)
    % Create default data
    
    portStruct = struct(...
        'polylineNumber',nan,...        % Number of polyline in polyline list
        'numberOfModes',1,...           % Number of modes to calculate
        'deembedDist',0,...             % Deembed distance, in simulation units
        'normImpedance',50,...          % Impedance to normalize to
        'power',1);                     % Power infused in simulation units
    
    portStruct = misc_validatePropStruct(portStruct,varargin);
    
end
