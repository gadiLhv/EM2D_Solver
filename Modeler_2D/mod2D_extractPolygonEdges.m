function edgeIdxs = mod2D_extractPolygonEdges(meshData,meshProps,polyToExtract)

    % Convert to polygon\node list, same as before.
    [polNodes,polEdges,~,~,~] = mod2D_polygonToFaceList({polyToExtract},[],meshProps.stitchingTolerance);
    
    % Extraction is done using using the .etri field
    edgePairs = meshData.etri;
    vert1 = meshData.vert(edgePairs(:,1),:);
    vert2 = meshData.vert(edgePairs(:,2),:);
    
%    % For each edge in the original polygon, match all the edge pairs
%    edgeVects = polNodes(polEdges(:,2),:) - polNodes(polEdges(:,1),:);
    
    [pointIdx1,whichEdge1,~,~] = mod2D_pointOnEdge(...
                                    vert1,...
                                    polNodes,...
                                    polEdges,...
                                    meshProps.stitchingTolerance);
                                    
    [pointIdx2,whichEdge2,~,~] = mod2D_pointOnEdge(...
                                    vert2,...
                                    polNodes,...
                                    polEdges,...
                                    meshProps.stitchingTolerance);
                                    
    % Detect points (equivalent to edge indexing) appearing in both
    
    [edgePoint1,~] = find(...
                    bsxfun(@eq,pointIdx1(:),pointIdx2(:).') ... 
                    & ...
                    bsxfun(@eq,whichEdge1(:),whichEdge2(:).'));
                    
    pointIdx1 = pointIdx1(edgePoint1);
    
    edgeIdxs = pointIdx1;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Debug: Draw polygon %
    %%%%%%%%%%%%%%%%%%%%%%%
%    figure;
%    
%    for edgeIdx = 1:size(polEdges,1)
%        plot(   polNodes(polEdges(edgeIdx,:),1),...
%                polNodes(polEdges(edgeIdx,:),2),...
%                '-b','linewidth',1);
%        hold on;
%    end
%
%    for edgeIdx = 1:numel(edgeIdxs);
%        cEdgePair = edgePairs(edgeIdxs(edgeIdx),:);
%        cVerts = meshData.vert(cEdgePair,:);
%        plot(cVerts(:,1),cVerts(:,2),'linewidth',2);
%    end
%
%    hold off;
%    
%    close(gcf);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End Debug: Draw polygon %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
end