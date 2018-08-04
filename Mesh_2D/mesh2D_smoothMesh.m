function smoothMesh = mesh2D_smoothMesh(meshToSmooth,meshProps)

          
% Smooth mesh
[ vert,...
  etri,...
  tria,...
  tnum] = ...
  smooth2(meshToSmooth.vert, ...
          meshToSmooth.etri, ...
          meshToSmooth.tria, ...
          meshToSmooth.tnum);

% Extract edge list and edge nodes for polylines
% 1. Determine how many polylines are there
lineIdxs = unique(meshToSmooth.lnum(:));

% 2. Iterate through assignements and create new lists
line = cell([numel(lineIdxs) 1]);
edge = [];
for lineIdx = 1:numel(lineIdxs)
  % Find all of the line assigments
  rowIdxs = find(meshToSmooth.lnum == lineIdxs(lineIdx));
  % From them, create a line assignment cell array
  line{lineIdx} = (1:numel(rowIdxs)).' + size(edge,1);
  % And create the edge list. 
  edge = [edge ; meshToSmooth.ltri(rowIdxs,:)];
  
  % All of these are referenced to the original vertex list
end



% Assign edges to polylines
[ ltri,...
  lnum] = ...
  mesh2D_assignEdgesToPolylines(...
        etri ,...
        vert,...
        line,...
        edge,...
        meshToSmooth.vert,...
        meshProps.stitchingTolerance)
          
% Update mesh structure
smoothMesh = meshToSmooth;
smoothMesh.vert = vert;
smoothMesh.etri = etri;
smoothMesh.tria = tria;
smoothMesh.tnum = tnum;
smoothMesh.ltri = ltri;
smoothMesh.lnum = lnum;
end