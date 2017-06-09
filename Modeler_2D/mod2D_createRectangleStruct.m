function rectStruct = mod2D_createRectangleStruct(r0,r1)
% Returns a polygon structure with "Rectangle" replaced 
% in the 'structType' field
% Inputs:
% r0 - (x,y) coordinates of start of rectangle
% r1 - (x,y) coordinates of end of retangle (opposite corner)
  
% Arrange minimum and maximum coor
x = [r0(1) r0(1) r1(1) r1(1)].';
y = [r0(2) r1(2) r1(2) r0(2)].';

% Check if start or end of square are the same
if((r0(1) == r1(1)) || (r0(2) == r1(2)))
  error('createRectangle',... % Error caller (id)
        'Start and end coordinates of rectangle cannot be the same. For a line, use ''createLine''');
end

rectStruct = mod2D_createPolygonStruct(x,y);
rectStruct.structType = 'Rectangle';
rectStruct.nParts = 1;

end