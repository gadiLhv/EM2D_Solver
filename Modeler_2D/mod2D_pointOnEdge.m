function [pointIdx,whichEdge,t,edgeLengths] = mod2D_pointOnEdge(vert,node,edge,minDistTH)
% Determine if the points in 'vert' are on any edge in 'node\edge'
% Inputs: 
% 1. vert - Vertices to snap to existing nodes. [Nverts,2]
% 2. node - Existing nodes. [Nnodes,2]
% 3. edge - Start and end node indexes. [Nedges,2]
% 4. minDistTH - minimum distance qualifier for vertex to be accepted as an existing node.
%
% Outputs:
% 1. pointIdx - Which of the vertices ("new" points) are on the edge
% 2. whichEdge - Which edge they are on
% 3. t - Normalized distance on the edge coordinate system (tangent direction)
% 4. edgeLengths - Lengths of edges 

% Take first and second nodes from each edege
node1 = node(edge(:,1),:);
node2 = node(edge(:,2),:);

% Determine distances between lines and individual nodes

% 1. Line difference vector
dr = node2 - node1;

% 2. Line coefficients according to line equation Ax + By + C = 0:
% dy*x-dx*y+dx*y1-dy*x1 = 0
A = dr(:,2);
B = -dr(:,1);
C = dr(:,1).*node1(:,2) - dr(:,2).*node1(:,1);

% 3. Determine distances between all vertices and all lines:
% |Ax + By + C|/sqrt(A^2+B^2)
d = bsxfun(@plus, bsxfun(@times,A.',vert(:,1)) + ...
                  bsxfun(@times,B.',vert(:,2)),...
                  C.');
d = bsxfun(@times,abs(d),1./sqrt((A.^2).' + (B.^2).'));

% Check normalized distance from node #1 
% 1. Determine edge lengths
L = sqrt(sum(dr.^2,2));

% 2. Normalize direction vector
dl = bsxfun(@times,dr,1./L);

% 3. Project (normalized coordinates, between 0 and 1)
node_offset = bsxfun(@minus,permute(vert,[1 3 2]),permute(node1,[3 1 2]));
t = sum(bsxfun(@times,node_offset,permute(dl,[3 1 2])),3);
t = bsxfun(@times,t,1./(L.'));

% Two conditions to snap a node to an edge:
% 1. First condition: distance smaller than threshold
distCond = d < minDistTH;


% 2. Second condition: Projection is withing bounds
t_de = bsxfun(@times,t,L.');
projCond = (t_de >= -minDistTH) & bsxfun(@le,t_de,L.' + minDistTH);
                
% Unify all three conditions
binSnapNodes = distCond & projCond;

% Extract only points who apply all 3 conditions
[pointIdx,whichEdge] = find(binSnapNodes);
% Also extract normalized distances.
t = t(pointIdx + (whichEdge-1)*size(t,1));

edgeLengths = L(whichEdge);

end