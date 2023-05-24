function l = cem2D_calcEdgeLengths(verts,triNodeIdxs)

% Extract all 3 nodes
r1 = verts(triNodeIdxs(:,1),:);
r2 = verts(triNodeIdxs(:,2),:);
r3 = verts(triNodeIdxs(:,3),:);

calcDist = @(p1,p2) sqrt((p2(:,1) - p1(:,1)).^2 + (p2(:,2) - p1(:,2)).^2);

l = [calcDist(r1,r2) calcDist(r3,r2) calcDist(r1,r3)];

end
