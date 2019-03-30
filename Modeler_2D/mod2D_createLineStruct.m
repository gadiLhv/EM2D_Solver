function lineStruct = mod2D_createLineStruct(r0,r1)
% Defines a single edge of a face as a port.
% Inputs:
% r0,r1 - 2 element vectors depicting the start and end point of the 

x = [r0(1) r1(1)].';
y = [r0(2) r1(2)].';

% Check if start or end of square are the same
if((r0(1) == r1(1)) && (r0(2) == r1(2)))
  error('createLine',... % Error caller (id)
        'Start and end coordinates of edge cannot be the same. For a edge, use ''createEdge''');
end
  
lineStruct.structType = 'Line';
lineStruct.x = x;
lineStruct.y = y;

% This field will contain the edges in the list.
lineStruct.edges = [];



end