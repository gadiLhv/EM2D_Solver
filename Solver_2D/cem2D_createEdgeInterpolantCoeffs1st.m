function ls = cem2D_createEdgeInterpolantCoeffs1st(vert,edge)

X = [vert(edge(:,1),1) vert(edge(:,2),1)];
Y = [vert(edge(:,1),2) vert(edge(:,2),2)];

ls = sqrt((X(:,2) - X(:,1)).^2 + (Y(:,2) - Y(:,1)).^2);

end