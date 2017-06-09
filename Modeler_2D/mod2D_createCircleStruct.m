function circStruct = mod2D_createCircleStruct(x0,y0,R)
% Returns a structure with parameters of a circle
% Fields:
% x0 - scalar, x coordinate of circle center
% y0 - scalar, y coordinate of circle center
% R - scalar, circle radius

circStruct = struct('x0',x0(1),'y0',y0(1),'R',R,'nParts',1,'structType','Circle');


end