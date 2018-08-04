function [ltri,lnum] = mesh2D_assignEdgesToPolylines(edgeVertIdxs,verts,linesAssignments,edges,node,minDistTH)

% Extract both vertices of face edges (etri)
vert1 = verts(edgeVertIdxs(:,1),:);
vert2 = verts(edgeVertIdxs(:,2),:);

%figure;
%for vIdx = 1:size(vert1,1)
%  plot([vert1(vIdx,1) vert2(vIdx,1)],[vert1(vIdx,2) vert2(vIdx,2)],'-b');
%  hold on;
%  plot(vert1(vIdx,1),vert1(vIdx,2),'.r','markersize',5);
%  plot(vert2(vIdx,1),vert2(vIdx,2),'.r','markersize',5);
%end
%hold off;

ltri = [];
lnum = [];

for lineIdx = 1:numel(linesAssignments)

  % Extract unmeshed edge list of the current polyline to be tested
  cLineEdgeIdxs = linesAssignments{lineIdx};
  cLineEdges = edges(cLineEdgeIdxs,:);
  
  % Check if both points 
  [vertIdx1,parentEdge1,~,~] = ...
    mod2D_pointOnEdge(  vert1,...         % Check if any of these
                        node,...          % Are on these edges
                        cLineEdges,...
                        minDistTH);
                        
  [vertIdx2,parentEdge2,~,~] = ...
    mod2D_pointOnEdge(  vert2,...         % Check if any of these
                        node,...          % Are on these edges
                        cLineEdges,...
                        minDistTH);

  % Check that both 1st and second vertex of each of these edges
  % are on the edge
  sameEdge = bsxfun(@eq,vertIdx1(:),vertIdx2(:).');
  sameParent = bsxfun(@eq,parentEdge1(:),parentEdge2(:).');
  
  sameEdgeAndParent = sameEdge & sameParent;
  [whichEdge1,whichEdge2] = find(sameEdgeAndParent);
  
%  figure;
%  plot(vert1(vertIdx1(whichEdge1),1),vert1(vertIdx1(whichEdge1),2),'ob');
%  hold on;
%  plot(vert2(vertIdx2(whichEdge2),1),vert2(vertIdx2(whichEdge2),2),'*r');
%  hold off;  
%  
  
  % Extract edge indices
  ltri = [ltri ; [edgeVertIdxs(vertIdx1(whichEdge1),1) edgeVertIdxs(vertIdx2(whichEdge2),2)]];
  % concatenate them into line 
  lnum = [lnum ; lineIdx*ones([numel(whichEdge1) 1])];
  
end

end