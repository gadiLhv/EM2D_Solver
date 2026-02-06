function [edgeIdxs,nEdges,edge2vert] = mesh2D_createEdgeIndexing(meshData)

vert_m = meshData.vert;
tris = meshData.tria;
% Convert vertex list to edge list
edge_m = [[tris(:,2) tris(:,1)] ; [tris(:,3) tris(:,2)] ; [tris(:,1) tris(:,3)]];

% Determine unique edges
edge_sorted = sort(edge_m,2);
[uEdges,uIdxs,eIdxs] = unique(edge_sorted,'rows');
% uEdges - List of unique edge pairs
% uIdxs - A list of the unique edges extracted from the entire edge list
% eIdxs - Assignment of edges to triangles. Namely, edge triplets.

nEdges = size(uEdges,1);

% Now I have a list that tells the triangle edge indexing. Namely, edge
% 'e_t' in triangle 't' is edge 'e' in the global edge list

% However, the list was 3 columns concatenated. Let's bring that back
edgeIdxs = reshape(eIdxs,[],3);

edge2vert = edge_m(uIdxs,:);

end

