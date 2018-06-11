function edgeCouples = cem2D_findEdgeNodes(faceTriplets,allEdges)
  
% Recover all edge nodes corresponding to the
% extracted triangles:
% 1. Determine if either of the edge nodes is present in a triangle
edgeNode1map = bsxfun(@eq,allEdges(:,1),permute(faceTriplets,[3 1 2]));
edgeNode2map = bsxfun(@eq,allEdges(:,2),permute(faceTriplets,[3 1 2]));
             
% 2. Determine if both of them are present in a triangle
areBothNodesIn = max(edgeNode1map,[],3) & max(edgeNode2map,[],3);

% 3. Mask node mappings
edgeNodeMap = bsxfun(@and,areBothNodesIn,edgeNode1map) | ...
              bsxfun(@and,areBothNodesIn,edgeNode2map);
% Permute so vertex numbers will come first
edgeNodeMap = permute(edgeNodeMap,[3 2 1]);

% 3. Extract only nodes that consist an edge in this face.
whichIdx = find(edgeNodeMap);
[whichVert,whichTri,whichEdge] = ind2sub(size(edgeNodeMap),whichIdx);
whichEdge = reshape(whichEdge,2,[]).';
whichVert = reshape(whichVert,2,[]).';
whichTri = reshape(whichTri,2,[]).';

% 4. Determine vertices order
% Sort ascending
whichVert = sort(whichVert,2);
% determine 1-3 paiStart extracting rs and reverse order
binPairs13 = (whichVert(:,1) == 1) & (whichVert(:,2) == 3);
pairsToFlip = whichVert(binPairs13,:);
pairsToFlip = pairsToFlip(:,2:-1:1);
whichVert(binPairs13,:) = pairsToFlip;

% 5. Extract edge couples
whichIdx1 = sub2ind(size(faceTriplets),whichTri(:,1),whichVert(:,1));
whichIdx2 = sub2ind(size(faceTriplets),whichTri(:,2),whichVert(:,2));
edgeCouples = [faceTriplets(whichIdx1) ...
               faceTriplets(whichIdx2)];

%edgeCouples = allEdges(whichEdge(:,1),:);


end