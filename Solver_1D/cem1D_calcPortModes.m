function [U,eigVal,f_cutoff,eff_diel] = cem1D_calcPortModes(portStruct, meshData, materialList, materialAssignment, meshProps, simProps)
    
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
        segMaterial{numel(segOrder)} = segmentWeakestMaterial(...
                                        segPt1, ...
                                        segPt2, ...
                                        meshData, ...
                                        materialList, ...
                                        materialAssignment, ...
                                        meshProps);
        
        % Replace the current point with the opposite on this segment
        cPt = segPt2;
        
        % Define     point at infinity so the distance will never be minimal
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
            [U,k_xi] = cem1D_calcTEmodes(segVerts,meshData,meshProps,segMaterial,materialList,simProps);
        case 'TM'
            [U,k_xi] = cem1D_calcTMmodes(segVerts,meshData,meshProps,segMaterial,materialList,simProps);
    end
    
        % Assumes that mu_r is 1!! (for effective DIELECTIC constant)
    c0 = physical_constant('speed of light in vacuum');
    e0 = physical_constant('electric constant');
    m0 = physical_constant('mag. constant');
    f0 = units(simProps.freqUnits,'Hz',simProps.fSim);
    k0 = 2*pi*f0/c0;
    
    eigVal = k_xi;
    kz = sqrt(k0^2 - k_xi.^2);
    f_cutoff = c0*k_xi/(2*pi);
    er_eff = c0*kz./(2*pi*f0)
   
end

function segMaterial = segmentWeakestMaterial(p1,p2,meshData, materialList, materialAssignment, meshProps);
    % Assign automatically with strongest materials, for start;
    segMaterial = 'PEC';
    segMaterialProps = cem2D_createMaterialDefs('name','PEC','type','PEC');
    minDistTH = meshProps.stitchingTolerance;
    % Iterate through face edges. First one is always bounding box
    for faceIdx = 1:numel(meshData.face)
        cEdges = meshData.edge(meshData.face{faceIdx},:);
        nodes1 = meshData.node(cEdges(:,1),:);
        nodes2 = meshData.node(cEdges(:,2),:);
        
        [pointIdx,whichEdge,~,~] = mod2D_pointOnEdge([p1 ; p2],...
                            meshData.node,...
                            meshData.edge(meshData.face{faceIdx},:),...     
                            minDistTH);
                                                
        % Assign whichEdge to the larger llist
        allEdgeIdxs = meshData.face{faceIdx};
        whichEdge = allEdgeIdxs(whichEdge);
        %%%%%%%%%
        % DEBUG %
        %%%%%%%%%
%        figure;
%        for cEdgeIdx = meshData.face{faceIdx}.';
%            cEdge = meshData.edge(cEdgeIdx,:);
%            cPts = meshData.node(cEdge.',:);
%            
%            plot(cPts(:,1),cPts(:,2),'-b');
%            hold on;
%        end
%        for cEdgeIdx = whichEdge.'
%            cEdge = meshData.edge(cEdgeIdx,:);
%            cPts = meshData.node(cEdge.',:);
%            
%            plot(cPts(:,1),cPts(:,2),'-r','linewidth',2);
%        end
%        plot(p1(1),p1(2),'.g','markersize',15);
%        plot(p2(1),p2(2),'.g','markersize',15);
%        hold off;
%        close(gcf);
        %%%%%%%%%
        % DEBUG %
        %%%%%%%%%
        
        % The two points should appear on the same edge twice
        maxPtsOnEdge = max(sum(bsxfun(@eq,whichEdge,whichEdge.'),2));
        % If the point list is empty, check that both points are on the shape
        if (~isempty(pointIdx)) && (maxPtsOnEdge == 2)
            cMaterialName = materialAssignment{faceIdx};
            cMaterialProps = cem2D_getMaterialPropsFromName(cMaterialName,materialList);
            
            [segMaterial,segMaterialProps] = chooseWeakestMaterial(segMaterial,segMaterialProps,cMaterialName,cMaterialProps);
        end
    end
end

function [materialName,materialProps] = chooseWeakestMaterial(material1name,material1defs,material2name,material2defs)
    materialRating = {'PEC','PMC','metal','normal'};
    get1Rating = @(matType) strcmp(material1defs.type,matType);
    get2Rating = @(matType) strcmp(material2defs.type,matType);
    
    rating1 = find(cellfun(get1Rating,materialRating));
    rating2 = find(cellfun(get2Rating,materialRating));
    
    if rating1 < rating2
        materialName = material2name;
        materialProps = material2defs;
    else
        materialName = material1name;
        materialProps = material1defs;
    end
end
