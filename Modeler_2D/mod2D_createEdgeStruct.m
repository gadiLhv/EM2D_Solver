function edgeStruct = mod2D_createEdgeStruct(r0,r1)
% Defines a single edge of a face as a port.
% Inputs:
% r0,r1 - 2 element vectors depicting the start and end point of the 

x = [r0(1) r1(1)].';
y = [r0(2) r1(2)].';

% Check if start or end of square are the same
if((r0(1) == r1(1)) && (r0(2) == r1(2)))
  error('createEdge',... % Error caller (id)
        'Start and end coordinates of edge cannot be the same. For a edge, use ''createEdge''');
end
  
edgeStruct = mod2D_createPolygonStruct(x,y);
edgeStruct.structType = 'Edge';

% In the beginning the edge has only one part.

% This field will contain the edges in the list.
edgeStruct.edges = [];



end