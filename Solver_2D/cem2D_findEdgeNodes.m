function edgeCouples = cem2D_findEdgeNodes(faceTriplets,allEdges)
  
% Recover all edge nodes corresponding to the
% extracted triangles:
% 1. Determine if either of the edge nodes is present in a triangle
edgeNode1map = bsxfun(@eq,allEdges(:,1),permute(faceTriplets,[3 1 2]));
edgeNode2map = bsxfun(@eq,allEdges(:,2),permute(faceTriplets,[3 1 2]));

% 2. Determine if both of them are present in a triangle
areBothNodesIn = max(edgeNode1map,[],3) & max(edgeNode2map,[],3);

% 3. Extract only nodes that consist an edge in this face.
[edgesInFace,~] = find(areBothNodesIn);
edgeCouples = allEdges(edgesInFace,:);

end