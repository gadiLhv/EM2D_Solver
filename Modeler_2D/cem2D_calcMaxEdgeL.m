function maxEdgeL = cem2D_calcMaxEdgeL(bBox)
% maxEdgeL = cem2D_calcMaxEdgeL(bBox)
%
% This function receives the bounding box of the entire structure and returns
% it's diagonal. This value is then sent to 'lfshfn2'. The resulting mesh then
% has new nodes only on the edges of the faces, rather than new nodes inside 
% the patch.

x1 = bBox.x(1);
x2 = bBox.x(3);

y1 = bBox.y(1);
y2 = bBox.y(3);

maxEdgeL = sqrt((x1-x2).^2 + (y1-y2).^2);

end