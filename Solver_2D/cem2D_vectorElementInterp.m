function fI = cem2D_vectorElementInterp(vert,triTriplets,edgeTriplets,edgeVals)
    % Creates one-time interpolant using a set of points and given values at
    % these points.

    dTri = delaunay(X(:),Y(:));

    Xv = [vert(triTriplets(:,1),1) ;

    % Calculate the coefficients for the scalar interpolant
    [a,b,c,Det] = cem2D_createInterpolantCoeffs1st(vert,triTriplets)

    % Calculate edge lengths
    eL = cem_calcEdgeLengths(vert,triTriplets);

    % Define interpolant
    fI = @(xI,yI) triInterpolant_1st(X,Y,f,dTri,a,b,c,Det,xI,yI);

end

function fI = vectInterpolant_1st(vert,dTri,a,b,c,Det,xI,yI);
    vert_m = [r1 ; r2 ; r3];
    [a,b,c,Det] = cem2D_createInterpolantCoeffs1st(vert_m,[1 2 3]);
    l1 = sqrt(sum((r2 - r1).^2,2));
    l2 = sqrt(sum((r3 - r2).^2,2));
    l3 = sqrt(sum((r1 - r3).^2,2));

    Li = @(i,x,y) (1/(2*Det))*(a(i) + b(i)*x + c(i)*y);
    dLidx = @(i,x,y) b(i);
    dLidy = @(i,x,y) c(i);

    N1x = @(x,y) l1*(Li(1,x,y)*b(2) - Li(2,x,y)*b(1));
    N1y = @(x,y) l1*(Li(1,x,y)*c(2) - Li(2,x,y)*c(1));

    N2x = @(x,y) l2*(Li(2,x,y)*b(3) - Li(3,x,y)*b(2));
    N2y = @(x,y) l2*(Li(2,x,y)*c(3) - Li(3,x,y)*c(2));

    N3x = @(x,y) l3*(Li(3,x,y)*b(1) - Li(1,x,y)*b(3));
    N3y = @(x,y) l3*(Li(3,x,y)*c(1) - Li(1,x,y)*c(3));



end


function fI = triInterpolant_1st(vert,f,dTri,a,b,c,Det,xI,yI)

    % Initialize
    fI = zeros(size(xI));

    % Find in which triangle is each
    triIdx = cem2D_pointInTri(vert,dTri,xI,yI);

    % Put NaNs for points not in mesh (perhaps zero is a better choice)
    fI(isnan(triIdx)) = NaN;

    % Note all points in mesh
    inMeshIdxs = find(~isnan(triIdx));

    % Construct interpolants for the points inside the mesh.

    % Extract ABCs and Determinant
    ai = a(triIdx(inMeshIdxs),:);
    bi = b(triIdx(inMeshIdxs),:);
    ci = c(triIdx(inMeshIdxs),:);
    Di = Det(triIdx(inMeshIdxs));

    % Extact weights
    triTriplets = dTri(triIdx(inMeshIdxs),:);
    wi = [f(triTriplets(:,1)) f(triTriplets(:,2)) f(triTriplets(:,3))];

    % Calculate weighted sum and place in result vector
    fI(inMeshIdxs) = (1./Di).*sum(0.5*(ai + xI(inMeshIdxs).*bi + yI(inMeshIdxs).*ci).*wi,2);

end

