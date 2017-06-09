function [node,connectivity,faceAssignment] = mod2D_snapLooseNodes(node,connectivity,faceAssignment,minDistTH)

% Repeat until there are no more edges to snap


while(1)
  % Take first and second nodes from each edege
  node1 = node(connectivity(:,1),:);
  node2 = node(connectivity(:,2),:);

  % Determine distances between lines and individual nodes

  % 1. Line difference vector
  dr = node2 - node1;

  % 2. Line coefficients according to line equation Ax + By + C = 0:
  % dy*x-dx*y+dx*y1-dy*x1 = 0
  A = dr(:,2);
  B = -dr(:,1);
  C = dr(:,1).*node1(:,2) - dr(:,2).*node1(:,1);

  % 3. Determine distances between all nodes and all lines:
  % |Ax + By + C|/sqrt(A^2+B^2)
  d = bsxfun(@plus, bsxfun(@times,A.',node(:,1)) + ...
                    bsxfun(@times,B.',node(:,2)),...
                    C.');
  d = bsxfun(@times,abs(d),1./sqrt((A.^2).' + (B.^2).'));

  % Check distance from node_1
  % 1. Determine edge lengths
  L = sqrt(sum(dr.^2,2));

  % 2. Normalize direction vector
  dl = bsxfun(@times,dr,1./L);

  % 3. Project
  node_offset = bsxfun(@minus,permute(node,[1 3 2]),permute(node1,[3 1 2]));
  t = sum(bsxfun(@times,node_offset,permute(dl,[3 1 2])),3);
  t = bsxfun(@times,t,1./(L.'));

  % Three conditions to snap a node to an edge:
  % 1. First condition: distance smaller than threshold
  distCond = d < minDistTH;

  % 2. Second condition: Projection is withing bounds
  projCond = (t >= 0) & (t <= 1);

  % 3. Node is not one of the two that compose of the edge. 
  notInEdge  = ~( bsxfun(@eq,(1:size(node,1)).',connectivity(:,1).') | ...
                  bsxfun(@eq,(1:size(node,1)).',connectivity(:,2).'));
                  
  % Unify all three conditions
  binSnapNodes = distCond & projCond & notInEdge;


  [nodeToSnap,edgeToSnapTo] = find(binSnapNodes);
  
  
  % If, and only if, there are no more nodes to snap, then you are free to go
  if(isempty(nodeToSnap))
    break;
  end
  
  
  % It is unknown how many nodes need to be snapped to each
  % edge, hence this must be done sequentially (per edge)
  uniqueEdgeIdx = unique(edgeToSnapTo);

  % Iterate through edges
  cEdgeIdx = uniqueEdgeIdx(1);
  
  % Get current list of nodes to snap
  binEdge = edgeToSnapTo == cEdgeIdx;
  cNodeToSnapIdx = nodeToSnap(binEdge);
  
  % Sort nodes according to distance from first node
  ts = t(cNodeToSnapIdx + (cEdgeIdx-1)*size(t,1));
  [~,sortIdx] = sort(ts);
  % Create new set of nodes
  
  cNodeToSnapIdx = [connectivity(cEdgeIdx,1) ; ...
                    cNodeToSnapIdx(sortIdx) ; ... 
                    connectivity(cEdgeIdx,2)];
  
  % Connectivity list to concatenate
  cConn = [cNodeToSnapIdx(1:end-1) cNodeToSnapIdx(2:end)];
  
  % Face assignement list to concatenate
  cfaceAssignment = ones([size(cConn,1) 1])*faceAssignment(cEdgeIdx);
  
  % Update connectivity and face assignements
  connectivity = [connectivity(1:(cEdgeIdx-1),:) ; ...
                  cConn ; ...
                  connectivity((cEdgeIdx+1):end,:)];
  
  faceAssignment = [ faceAssignment(1:(cEdgeIdx-1),:) ; ...
                      cfaceAssignment ; ...
                      faceAssignment((cEdgeIdx+1):end,:)];
                      

end


end