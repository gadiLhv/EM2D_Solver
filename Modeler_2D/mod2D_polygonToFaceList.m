function [nodes,edges,faces,lines,fixesApplied] = mod2D_polygonToFaceList(polygons,polylines,uniqueDistTH)
    % [nodes,faces,edge] = polygonToFaceList(polygons)
    %  
    % edges - [Ne,3] edge indexing. Each edge is represented by it's connectivity.
    %         edge(:,1) - Edge indices. This is a legend for later use.
    %         edge(:,2:3) - First and second NODE indices. This, uniquely represents
    %         and edge.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Convert all polygons to node\connectivity lists and join to a node list %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allNode = [];
    allConnect = [];
    faceAssignments = [];
    
    nodeCounter = 0;
    
    %figure;
    %plot(0,0);
    %axHdl = gca;
    
    for polIdx = 1:numel(polygons)
        % Convert to node\connectivity list
        [x,y,cConn] = mod2D_polygon2connect(polygons{polIdx});
        
        cNode = [x y];
        
        % Concatenate to list
        allNode = [allNode ; cNode];
        
        % Attach currnet connectivity list
        allConnect = [allConnect ; (cConn + nodeCounter)];
        
        % Face assignements are simply all the current edges added
        faceAssignments = [ faceAssignments ; polIdx*ones([size(cConn,1) 1])];
        
        % Increment node counter
        nodeCounter = nodeCounter + size(cNode,1);
    end
    
    % Store the maximum number
    nFacesPolygons = max(faceAssignments(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Clip polylines with all polygons. Mock adds polylines as faces, %
    % for smoother integration.                                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for lineIdx = 1:numel(polylines)
        cLine = polylines{lineIdx};
        cLine = [cLine.x cLine.y];
        
        % Clip with all polygons before continuing. Skip bounding box.
        for polIdx = 2:numel(polygons)
            % Extract current polygon
            %    [x,y,cConn] = mod2D_polygon2connect(polygons{polIdx});
            [x,y] = polygon2patch(polygons{polIdx}.x,polygons{polIdx}.y);
            cPol = [x y];
            
            % Subtract op
            [subLine,~] = clipPolyline(cLine, cPol, 0);
            % Add Op
            [addLine,~] = clipPolyline(cLine, cPol, 1);
            
            % Add to single list
            allSegs = [[NaN NaN] ; addLine ; [NaN NaN] ; subLine ; [NaN NaN]];
            
            % Rebuild line by iteratively searching the first point
            cPoint = cLine(1,:); 
            clippedLine = [];
            
            % This is based on the assumption that the first and last point of the line
            % only appear once.
            while (sum(isnan(allSegs(:))) ~= numel(allSegs))
                % Mark segments
                segDelimeters = find(isnan(allSegs(:,1)));
                segDelimeters = [segDelimeters(1:(end-1)) segDelimeters(2:end)];
                
                
                % Find the row of this segment
                [whichIdx,~] = find(sum(bsxfun(@eq,cPoint,allSegs),2) == 2);
                
                % Find which segment is this
                whichSeg = find((whichIdx > segDelimeters(:,1)) & (whichIdx < segDelimeters(:,2)));
                
                % Determine if the point is the first or last in the segment, and arrange
                % accordingly.
                segStart = segDelimeters(whichSeg,1) + 1;
                segEnd = segDelimeters(whichSeg,2) - 1;
                
                isStart = sum(allSegs(segStart,:) == cPoint,2) == 2;
                if isStart
                    clippedLine = [clippedLine ; allSegs(segStart:segEnd,:)];
                    else
                        clippedLine = [clippedLine ; allSegs(segEnd:-1:segStart,:)];
                    end
                    
                    % Delete segment from the segment list
                    allSegs((segStart - 1):segEnd,:) = [];
                    
                    % Denote new point as current last one
                    cPoint = clippedLine(end,:);
                end
                
                
                % Remove duplicate lines
                [clippedLine,~,~] = mod2D_uniqueNodeByDist(clippedLine,uniqueDistTH);
                
                cLine = clippedLine;
                
            end
            
%            %%%%%%%%%%%%%%%%%%%%%%%%
%            % Debug: Show polygons %
%            % and current line     %
%            %%%%%%%%%%%%%%%%%%%%%%%%
%            figure('position',[852   163   560   420]);;
%            [x,y] = polygon2patch(polygons{1}.x,polygons{1}.y);
%            patch(x,y,'facecolor',[1 1 1],'edgecolor',[0 0 0]);
%            hold on;
%            for polIdx = 2:numel(polygons)
%                [x,y] = polygon2patch(polygons{polIdx}.x,polygons{polIdx}.y);
%                patch(x,y,'facecolor',rand([1 3]),'edgecolor',[0 0 0]);
%            end
%            plot(cLine(:,1),cLine(:,2),'-r','linewidth',3);
%            plot(cLine(:,1),cLine(:,2),'og','linewidth',3);
%            hold off;
%            close(gcf);
            
            % Count node indexes
            cNodeIdxs = (size(allNode,1) + (1:size(cLine,1))).';
            
            % Insert nodes and edges
            cEdges = [cNodeIdxs(1:(end-1)) cNodeIdxs(2:end)];
            
            % Add to list
            allNode = [allNode ; cLine];
            allConnect = [allConnect ; cEdges];
            
            % Effectively add as face, for node snapper function. Will be removed later.
            faceAssignments = [faceAssignments ; (nFacesPolygons + lineIdx)*ones([size(cEdges,1) 1])];
        end
        
        cNumNodes = size(allNode,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Remove duplicate nodes %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        [ allNode,allConnect,faceAssignments] = mod2D_removeDuplicateNodes(...
        allNode,...            % All nodes
        allConnect,...         % Current connectivity list
        faceAssignments,...    % Current face assignement
        uniqueDistTH);         % Distance to consider the same node
        
        fixesApplied{1} = sprintf('%d nodes removed as duplicate. Dth = %.3e\n',cNumNodes-size(allNode,1),uniqueDistTH);
        
        cNumEdges = size(allConnect,1);
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Snapping loose nodes %
        %%%%%%%%%%%%%%%%%%%%%%%%
        % In case (again, due to fine thresholds) that nodes are "off edge", snap new
        % edges.
        [allNode,allConnect,faceAssignments] = mod2D_snapLooseNodes(...
        allNode,...
        allConnect,...
        faceAssignments,...
        uniqueDistTH);
        
        fixesApplied{2} = sprintf('%d nodes snapped to close edges. Dth = %.3e\n',size(allConnect,1)-cNumEdges+1,uniqueDistTH);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Re-establish face and line assignements %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Entire possible edge list
        edgeIdx = (1:size(allConnect,1)).';
        % Assign to faces
        faces = cell([numel(polygons) 1]);
        for faceIdx = 1:numel(polygons)
            faces{faceIdx} = edgeIdx(faceAssignments == faceIdx);
        end
        
        % Assign to lines
        lines = cell([numel(polylines) 1]);
        for lineIdx = (numel(polygons) + (1:numel(polylines)))
            lines{lineIdx - numel(polygons)} = edgeIdx(faceAssignments == lineIdx);
        end
        
        % And completely remove the polylines from the list
        faceAssignments(faceAssignments > nFacesPolygons) = [];
        
        % Store to output variables
        nodes = allNode;
        %connectivity = allConnect;
        
        % Each unique pair in the 'allConnect' list represents an edge. 
        % For now assume it's each pair....
        edges = allConnect;
        
        
        % Face assignment is per-edge
        
        
    end
