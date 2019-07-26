function segMaterial = misc_segmentWeakestMaterial(p1,p2,meshData, materialList, materialAssignment, meshProps);
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
            
            [segMaterial,segMaterialProps] = misc_chooseWeakestMaterial(segMaterial,segMaterialProps,cMaterialName,cMaterialProps);
        end
    end
end
