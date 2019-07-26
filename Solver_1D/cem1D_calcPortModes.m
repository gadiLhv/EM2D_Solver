function portStruct = cem1D_calcPortModes(portStruct, meshData, materialList, materialAssignment, meshProps, simProps)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                              %%
    %% Assign materials to segments %%
    %%                              %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Extract line and vertices.
    lineNum = portStruct.polylineNumber;
    lineVertPairs = meshData.ltri(meshData.lnum == lineNum,:);
    lineR1 = meshData.vert(lineVertPairs(:,1),:);
    lineR2 = meshData.vert(lineVertPairs(:,2),:);
        
    allPoints = [lineR1 ; lineR2];
    
    % Detect direction of line
    % 1. Find the farthermost points
    % Create distance map
    distMap = sqrt(bsxfun(@minus,allPoints(:,1),allPoints(:,1).').^2 + ...
                   bsxfun(@minus,allPoints(:,2),allPoints(:,2).').^2);
    % Extract two maximum points. These are the edges
    [portLength,maxIdx] = max(distMap(:));
    
    % Give points meaning in larger list
    maxRow = mod(maxIdx - 1,size(allPoints,1)) + 1;
    maxCol = floor((maxIdx - 1)/size(allPoints,1)) + 1;
    % Choose arbitrarily the start and stop points
    p1 = allPoints(maxRow,:);
    p2 = allPoints(maxCol,:);
    segDir = p2 - p1;
    segNorm = ([[0 -1] ; [1 0]]*segDir(:)).';
    
    % Calculate port center for later
    portCenter = 0.5*(p1 + p2);
    % And normalize direction
    segNorm = segNorm/sqrt(sum(segNorm(:).^2,1));
    
    % Segment by segment, build adjacent material list. The material considered
    % as the WG structure is the "weakest" in the list. The order from strongest
    % to weakest:
    % 1. PEC
    % 2. PMC
    % 3. metal (finite conductivity\lossy)
    % 4. Dielectric.
    % Any contradiction (i.e. two dielectrics, two metals, etc) results in an
    % error.
    
    % 2. Now start arranging the segment piece by piece
    segOrder = [];
    segMaterial = {};
    cPt = p1;
    calcDist = @(x1,x2) sqrt((x1(:,1) - x2(:,1)).^2 + (x1(:,2) - x2(:,2)).^2);
    fLineR1 = lineR1;
    fLineR2 = lineR2;
    for pIdx = 1:size(lineR1,1)
        % Create distance list
        [minDists,closestSeg] = min([calcDist(cPt,fLineR1) calcDist(cPt,fLineR2)]);
        
        [~,whichIsMin] = min(minDists(:));
        % Write down this current segment
        closestSeg = closestSeg(whichIsMin);
        segOrder = [segOrder ; closestSeg];
        
        segPt1 = cPt;
        segPt2 = fLineR1(closestSeg,:)*(whichIsMin == 2) + fLineR2(closestSeg,:)*(whichIsMin == 1);
        % Take the two segment points and look for all overlapping materials.
        segMaterial{numel(segOrder)} = misc_segmentWeakestMaterial(...
                                        segPt1, ...
                                        segPt2, ...
                                        meshData, ...
                                        materialList, ...
                                        materialAssignment, ...
                                        meshProps);
        
        % Replace the current point with the opposite on this segment
        cPt = segPt2;
        
        % Define point at infinity so the distance will never be minimal
        fLineR1(closestSeg,:) = inf;
        fLineR2(closestSeg,:) = inf;
        
    end
    
    segVerts = lineVertPairs(segOrder,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Debug: Test segment ordering %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    figure;
%    for segIdx = 1:numel(segOrder)
%        cSeg = segOrder(segIdx);
%        p1 = lineR1(cSeg,:);
%        p2 = lineR2(cSeg,:);
%        
%        pMean = 0.5*(p1 + p2);
%        
%        plot([p1(1) ; p2(1)],[p1(2) ; p2(2)],'-b');
%        hold on;
%        plot([p1(1) ; p2(1)],[p1(2) ; p2(2)],'or');
%        textHdl = text(pMean(1),pMean(2),[num2str(segIdx) ':' num2str(cSeg) ' - ' segMaterial{segIdx}]);
%        set(textHdl,'fontsize',14);
%    end
%    hold off;
%    set(gca,'fontsize',14);
%    xlabel('x [mm]','fontsize',16);
%    ylabel('y [mm]','fontsize',16);
%    
%    close(gcf);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End Debug: Test segment ordering %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                         %%
    %% Build Eigenvalue Matrix %%
    %%                         %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%    
    switch simProps.polarizationType
        case 'TE'
            [U,k_xi,segVerts] = cem1D_calcTEmodes(segVerts,meshData,meshProps,segMaterial,materialList,simProps);
        case 'TM'
            [U,k_xi,segVerts] = cem1D_calcTMmodes(segVerts,meshData,meshProps,segMaterial,materialList,simProps);
    end
    
        % Assumes that mu_r is 1!! (for effective DIELECTIC constant)
    c0 = physical_constant('speed of light in vacuum');
    e0 = physical_constant('electric constant');
    m0 = physical_constant('mag. constant');
    f0 = units(simProps.freqUnits,'Hz',simProps.fSim);
    k0 = 2*pi*f0/c0;
    
    eigVal = k_xi;
    portStruct.f_cutoff = units('Hz',simProps.freqUnits,c0*k_xi/(2*pi));
    portStruct.portSegments = segVerts;
    portStruct.portModes = U;
    portStruct.portDirection = segNorm;
    portStruct.portCenter = portCenter;
    portStruct.segMaterial = segMaterial;
end
