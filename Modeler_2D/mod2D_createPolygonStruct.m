function polStruct = mod2D_createPolygonStruct(X,Y)
% Returns a structure with the (x,y) coordinates of the 
% Fields:
% x - X coordinates, aligned in a single column
% y - Y coordinates, aligned in a single column
% nParts - Number of closed loops that composes polygons
% structType - string. Name of the element type
%
% No need to repeat first coordinate twice

polStruct = struct('x',X(:),'y',Y(:),'nParts',1,'structType','Polygon');

% TODO: Determine number of loops, and write them in structure



end