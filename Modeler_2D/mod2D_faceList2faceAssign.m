function [connectivity,faceAssignment] = mod2D_faceList2faceAssign(edge,face)

faceAssignment = [];
connectivity = [];
% Concatenate face assignement
for faceIdx = 1:numel(face)
  cFace = face{faceIdx};
  
  % Current face assignment is simply the index of the current face
  faceAssignment = [faceAssignment ; ones([size(cFace,1) 1])*faceIdx];
  
  % Current connectivity is the extraction of all the relevant edges
  connectivity = [connectivity ; edge(cFace,:)];
end


end