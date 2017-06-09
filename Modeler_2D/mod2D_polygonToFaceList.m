function [nodes,edges,faces,fixesApplied] = mod2D_polygonToFaceList(polygons,uniqueDistTH)
% [nodes,faces,edge] = polygonToFaceList(polygons)
%  
% edges - [Ne,3] edge indexing. Each edge is represented by it's connectivity.
%         edge(:,1) - Edge indices. This is a legend for later use.
%         edge(:,2:3) - First and second NODE indices. This, uniquely represents
%         and edge.
  
% Convert all polygons to node\connectivity lists and join to a node list
allNode = [];
allConnect = [];
faceAssignments = [];

nodeCounter = 0;
for polIdx = 1:numel(polygons)
  % Convert to node\connectivity list
  [x,y,cConn] = mod2D_polygon2connect(polygons{polIdx});
  
  cNode = [x y];
  
  % Concatenate to list
  allNode = [allNode ; cNode];
  
  % Attach currnet connectivity list
  allConnect = [allConnect ; (cConn+nodeCounter)];
  
  % Face assignements are simply all the current edges added
  faceAssignments = [ faceAssignments ; polIdx*ones([size(cConn,1) 1])];
  
  % Increment node counter
  nodeCounter = nodeCounter + size(cNode,1);
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

fprintf('%d nodes removed as duplicate. Dth = %.3e\n',cNumNodes-size(allNode,1),uniqueDistTH);
                        
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

fprintf('%d nodes snapped to close edges. Dth = %.3e\n',size(allConnect,1)-cNumEdges+1,uniqueDistTH);
% Re-establish face assignements
edgeIdx = (1:size(allConnect,1)).';
% Now assign to faces
faces = cell([numel(polygons) 1]);
for faceIdx = 1:numel(polygons)
  faces{faceIdx} = edgeIdx(faceAssignments == faceIdx);
end


% Store to output variables
nodes = allNode;
%connectivity = allConnect;

% Each unique pair in the 'allConnect' list represents an edge. 
% For now assume it's each pair....
edges = allConnect;


% Face assignment is per-edge


end