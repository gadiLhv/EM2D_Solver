function [K,b] = cem2D_createKmatBvect_WG_Port(meshData,portStruct,materialList,simProps,meshProps,f_sim)
    
    % Vertices in which the field eigenvectors 
    segVerts = portStruct.portSegments;
    % Cutoff frequency in simulation units
    fc = portStruct.f_cutoff;
    fc = units(simProps.freqUnits,'Hz',fc);
    % Port simulation frequency
    f0 = units(simProps.freqUnits,'Hz',f_sim);
    % Constant
    c0 = physical_constant('speed of light in vacuum');
    
    % Calculate the propagation wave number
    kz = (2*pi/c0)*sqrt(f0^2 - fc.^2);
    
    % De-embedding distance
    deDist = units(simProps.lengthUnits,'m',portStruct.deembedDist);
    
    Ncells = size(meshData.vert,1);
    
    % Initialize solution matrix\vector. Later on this should be a sparse matrix
    K = zeros([Ncells*[1 1] portStruct.numberOfModes]);
    b = zeros([Ncells 1 portStruct.numberOfModes]);
    
    % Direction (for boundary conditions, later)
    nx = portStruct.portDirection(1)*portStruct.propDirection;
    ny = portStruct.portDirection(2)*portStruct.propDirection;
    
    % Material for each segment
    segMaterial = portStruct.segMaterial;
    
    % Calculate segment lengths
    pt1 = units(simProps.lengthUnits,'m',meshData.vert(segVerts(:,1),:));
    pt2 = units(simProps.lengthUnits,'m',meshData.vert(segVerts(:,2),:));
    l = sqrt(sum(abs(pt2 - pt1).^2,2));
    
    % Maps segments by 
    [uniqueSegVerts,~,segVertsMapping] = unique(segVerts(:));
    segVertsMapping = reshape(segVertsMapping,[],2);
    
    
    for modeIdx = 1:portStruct.numberOfModes
        % Extract quantities one by one
        gamma = -1i*kz(modeIdx);
        q = -2*1i*kz(modeIdx)*portStruct.portModes(:,:,modeIdx)*exp(-1i*kz(modeIdx)*deDist);
        
        cK = zeros(numel(uniqueSegVerts)*[1 1]);
        cb = zeros([numel(uniqueSegVerts) 1]);
        for segIdx = 1:size(segVerts,1)
            % Current nodes\vertices to deal with in the sub matrix
            segRange = segVertsMapping(segIdx,:);
            
            % Add this sub matrix to the port matrix
            cK(segRange,segRange) = cK(segRange,segRange) + gamma*l(segIdx)*[1 2 ; 2 1]/6;
            cb(segRange) = cb(segRange) + q(segIdx,:).'*l(segIdx)/2;
            
        end
        
        K(uniqueSegVerts,uniqueSegVerts,modeIdx) = cK;
        b(uniqueSegVerts,1,modeIdx) = cb;
    end

end