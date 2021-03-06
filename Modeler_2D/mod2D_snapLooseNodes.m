function [node,connectivity,faceAssignment] = mod2D_snapLooseNodes(node,connectivity,faceAssignment,minDistTH)

% Repeat until there are no more edges to snap


while(1)
    
  [nodeToSnap,edgeToSnapTo,t,edgeLengths] = ...
    mod2D_pointOnEdge(  node,...
                        node,...
                        connectivity,...
                        minDistTH);
  
  
  % Exclude "self" nodes (t == 0 or t==1), with respect to threshold
  selfNode = (abs(t).*edgeLengths <= minDistTH) | (abs(t - 1).*edgeLengths <= minDistTH);
  nodeToSnap(selfNode) = [];
  edgeToSnapTo(selfNode) = [];
  t(selfNode) = [];
  
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
  ts = t(binEdge);
  
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