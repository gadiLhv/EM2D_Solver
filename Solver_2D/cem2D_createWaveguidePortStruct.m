function portStruct = cem2D_createWaveguidePortStruct(varargin)
    % Create default data
    
    portStruct = struct(...
        'polylineNumber',nan,...        % Number of polyline in polyline list
        'numberOfModes',1,...           % Number of modes to calculate
        'segVerts',[],...               % Segment vertex pairs in given mesh
        'portModes',[],...              % THE port modes. Given in default, V/m or A/m
        'f_cutoff',[],...               % Cutoff frequency given in simulation units
        'deembedDist',0,...             % Deembed distance, in simulation units
        'normImpedance',50,...          % Impedance to normalize to
        'power',1);                     % Power infused in simulation units
    
    portStruct = misc_validatePropStruct(portStruct,varargin);
    
end
