function [nodes,edges,faces,edgeList,fixesApplied] = mod2D_addEdgesToFaceList(nodes,edges,faces,edgeList,minDistTH)


% Assign node numbers for new nodes
newNodeIdxs = size(nodes,1) + (1:numel(edgeList)*2)).';
% Assign edges of new list
newEdges = reshape(newNodeIdxs,[],2).';

newEdgeIdxs = size(edges,1) + (1:size(newEdges,1));

newEdgeNodes = [];
% Edge by edge, check if points are on edge
for edgeIdx = 1:numel(edgeList)
  cEdge = edgeList{edgeIdx};
 
  % Concatenate these nodes to the larger node list
  newEdgeNodes = [newEdgeNodes ; [cEdge.x cEdge.y]]; 
  
end

%% Build full node list.
%nodes = [nodes ; edgeNodes];
%allNodeIdxs = (1:size(nodes,1)).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove duplicate nodes, and if necessary, remove entire "new" edge %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[t_new,t_exist,L_new,L_exist] = mod2D_findEdgeIntersections( ...
                                  edgeNodes(1:2:end,:), ...
                                  edgeNodes(2:2:end,:), ...
                                  nodes(edges(:,1),:), ...
                                  nodes(edges(:,2),:));

% First, find which points ARE actually on the existing edges
isNew1stPoint = abs(bsxfun(@times,t_new,L_new) <= minDistTH);
isNew2ndPoint = abs(bsxfun(@times,1 - t_new,L_new) <= minDistTH);

% Only check for first point in edges, as all loops need to close.
isExist1stPoint = abs(bsxfun(@times,t_exist,L_exist) <= minDistTH);
% Check second edge to see if entire edge need to be replaced.
isExist2ndPoint = abs(bsxfun(@times,1 - t_exist,L_exist) <= minDistTH);

% Two cases:
% 1. Only one point need to be replaced.
replace1stPointBy1stPoint = isNew1stPoint & isExist1stPoint;
replace1stPointBy2ndPoint = isNew1stPoint & isExist2ndPoint;
replace2ndPointBy1stPoint = isNew2ndPoint & isExist1stPoint;
replace2ndPointBy2ndPoint = isNew2ndPoint & isExist2ndPoint;

% 2. Replace entire edge.
replaceEdge12 = replace1stPointBy1stPoint & replace2ndPointBy2ndPoint;
replaceEdge21 = replace1stPointBy2ndPoint & replace2ndPointBy1stPoint;


% The two are mutually exclusive. Namely, an entire edge excludes 
% a single point replacement.
replace1stPoint = replacePoint & ~replaceEdge;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replace one point in edge %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Find nodes to replace
[whichNewIdx,whichExistIdx] = find(replacePoint);

% In case node is common to several edges, make sure that the node to be
% replaced is unique.
[whichNewIdx,orig2unique,~] = unique(whichNewIdx);
whichExistIdx = whichExistIdx(orig2unique);

% And assign in the new edges list
newEdges(whichNewIdx) = edges(whichExistIdx,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check intersections with existing edges %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                                  

% Is on edge
doesIntersect = (t_new > 0) & (t_new < 1) & (t_exist > 0) & (t_exist < 1);

% Is an existing node
isFirstPoint = abs(t.*edgesLength) <= minDistTH;
isSecondPoint = abs((1 - t).*edgesLength) <= minDistTH;

% The two are mutually exclusive. Meaning, that if
% A point is an existing node, write it down and replace.
doesIntersect(isFirstPoint | isSecondPoint) = false(1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the edge points are the same as exiting nodes, simply replace the %
%    corrdinates and denote it.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change existing points
newEdges(isFirstPoint,1) = edges(whichEdge,1);
newEdges(isSecondPoint,2) = edges(whichEdge,2);

% Rebuild the node list

% 2. Remove duplicates
[uNode,uIdx,u2oIdx] = mod2D_uniqueNodeByDist(node,distTH)

% 3. Reconstruct index list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the happy edge received an additional node, split it to two edges %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end