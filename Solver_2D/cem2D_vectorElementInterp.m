function fI = cem2D_vectorElementInterp(vert,triTriplets,edgeTriplets,edgeVals)
    % Creates one-time interpolant using a set of points and given values at
    % these points.

    % Calculate the coefficients for the scalar interpolant
    [a,b,c,Det] = cem2D_createInterpolantCoeffs1st(vert,triTriplets);

    % Calculate edge lengths
    eL = cem2D_calcEdgeLengths(vert,triTriplets);

    % Define interpolant
    fI = @(xI,yI) vectInterpolant_1st(vert,triTriplets,edgeTriplets,edgeVals,xI,yI);

end

function fI = vectInterpolant_1st(vert,dTri,eIdxs,Ee,xI,yI);
% vert - All nodes
% dTri - Triangle node assignments
% eIdxs - Edge triplet indexes
% f_e  - edge element coefficients
% xI,yI - Values to interpolate

    xI = xI(:);
    yI = yI(:);

    % Initialize
    fI = ones([numel(xI) 2])*NaN;

    % Find in which triangle is each
    triIdx = cem2D_pointInTri(vert,dTri,xI,yI);

    % Put NaNs for points not in mesh (perhaps zero is a better choice)
    fI(isnan(triIdx),:) = NaN;

    % Note all points in mesh
    inMeshIdxs = find(~isnan(triIdx));

    % Interpolation coefficients
    [a,b,c,Det] = cem2D_createInterpolantCoeffs1st(vert,dTri);

    % lengths of all edges
    eL = cem2D_calcEdgeLengths(vert,dTri);

    % Extract ABCs and Determinant
    ai = a(triIdx(inMeshIdxs),:);
    bi = b(triIdx(inMeshIdxs),:);
    ci = c(triIdx(inMeshIdxs),:);
    Di = Det(triIdx(inMeshIdxs));

    l1 = eL(triIdx(inMeshIdxs),1);
    l2 = eL(triIdx(inMeshIdxs),2);
    l3 = eL(triIdx(inMeshIdxs),3);

    xI = xI(inMeshIdxs);
    yI = yI(inMeshIdxs);

    Li = @(i,x,y) (ai(:,i) + bi(:,i).*x + ci(:,i).*y)./(2*Di);

    N1x = l1.*(Li(1,xI,yI).*bi(:,2) - Li(2,xI,yI).*bi(:,1))./(2*Di);
    N1y = l1.*(Li(1,xI,yI).*ci(:,2) - Li(2,xI,yI).*ci(:,1))./(2*Di);

    N2x = l2.*(Li(2,xI,yI).*bi(:,3) - Li(3,xI,yI).*bi(:,2))./(2*Di);
    N2y = l2.*(Li(2,xI,yI).*ci(:,3) - Li(3,xI,yI).*ci(:,2))./(2*Di);

    N3x = l3.*(Li(3,xI,yI).*bi(:,1) - Li(1,xI,yI).*bi(:,3))./(2*Di);
    N3y = l3.*(Li(3,xI,yI).*ci(:,1) - Li(1,xI,yI).*ci(:,3))./(2*Di);


    Ee_1 = Ee(eIdxs(triIdx(inMeshIdxs),1));
    Ee_2 = Ee(eIdxs(triIdx(inMeshIdxs),2));
    Ee_3 = Ee(eIdxs(triIdx(inMeshIdxs),3));

    fxI = N1x.*Ee_1 + N2x.*Ee_2 + N3x.*Ee_3;
    fyI = N1y.*Ee_1 + N2y.*Ee_2 + N3y.*Ee_3;

    fI(inMeshIdxs,:) = [fxI fyI];
end

