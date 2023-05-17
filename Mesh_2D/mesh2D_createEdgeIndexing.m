function [edgeIdxs,nEdges] = mesh2D_createEdgeIndexing(meshData)

vert_m = meshData.vert;
tris = meshData.tria;
% Convert vertex list to edge list
edge_m = [[tris(:,2) tris(:,1)] ; [tris(:,3) tris(:,2)] ; [tris(:,1) tris(:,3)]];

edgeCent_x = [  0.5*(vert_m(tris(:,2),1) + vert_m(tris(:,1),1)) ; ...
                0.5*(vert_m(tris(:,3),1) + vert_m(tris(:,2),1)) ; ...
                0.5*(vert_m(tris(:,1),1) + vert_m(tris(:,3),1))];

edgeCent_y = [  0.5*(vert_m(tris(:,2),2) + vert_m(tris(:,1),2)) ; ...
                0.5*(vert_m(tris(:,3),2) + vert_m(tris(:,2),2)) ; ...
                0.5*(vert_m(tris(:,1),2) + vert_m(tris(:,3),2))];

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

% Edges have different directions in various triangles, hence it is required
% to store them again (maybe)
%edgeCent = [edgeCent_x(uIdxs) edgeCent_y(uIdxs)];

end

