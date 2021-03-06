function [a,b,c,Det] = cem2D_createInterpolantCoeffs1st(vert,tri)

X = [vert(tri(:,1),1) vert(tri(:,2),1) vert(tri(:,3),1)];
Y = [vert(tri(:,1),2) vert(tri(:,2),2) vert(tri(:,3),2)];


a = [ X(:,2).*Y(:,3) - Y(:,2).*X(:,3) ...
      X(:,3).*Y(:,1) - Y(:,3).*X(:,1) ...
      X(:,1).*Y(:,2) - Y(:,1).*X(:,2)];
b = [ Y(:,2) - Y(:,3) ...
      Y(:,3) - Y(:,1) ...
      Y(:,1) - Y(:,2)];
c = [ X(:,3) - X(:,2) ...
      X(:,1) - X(:,3) ...
      X(:,2) - X(:,1)];

Det = 0.5*(b(:,1).*c(:,2) - b(:,2).*c(:,1));

end