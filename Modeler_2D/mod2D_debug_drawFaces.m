function figHdl = mod2D_debug_drawFaces(node,edge,face,faceIdx)

% If face is a cell array, then this is after face assignment. If not, 
% assign faces first.
if ~iscell(face)
  % Re-establish face assignements
  edgeIdx = (1:size(edge,1)).';
  % Now assign to faces
  faces = cell([numel(unique(face)) 1]);
  for idx = 1:numel(faces)
    faces{idx} = edgeIdx(face == idx);
  end
  
  face = faces;
end



% Extract relevant edges
cEdges = edge(face{faceIdx},:);


cNodesX = [node(cEdges(:,1),1) node(cEdges(:,2),1)];
cNodesY = [node(cEdges(:,1),2) node(cEdges(:,2),2)];


plot(cNodesX,cNodesY,'.','markersize',25);
hold on;
for edgeIdx = 1:size(cEdges,1);
  
  
  arrowStop = node(cEdges(edgeIdx,2),:);
  arrowStart = node(cEdges(edgeIdx,1),:);
  arrowDir = arrowStop - arrowStart;
  
  plot([arrowStart(1) arrowStop(1)],[arrowStart(2) arrowStop(2)],'linewidth',3);
  quiver(arrowStart(:,1),arrowStart(:,2),arrowDir(:,1),arrowDir(:,2),'linewidth',3);
  


end

figHdl = gcf;
hold off;