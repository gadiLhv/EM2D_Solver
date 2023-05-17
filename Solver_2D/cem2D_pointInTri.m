function triIdx = cem2D_pointInTri(vert,dTri,Xp,Yp)
    Xv = vert(:,1);
    Yv = vert(:,2);

    % Define the maximum (approximate) memory usage for this function
    maxChunkSizeMB = 1000;

    % Determine in which trinagle of the triangle is each coordinate in
    R1 = [Xv(dTri(:,1)) Yv(dTri(:,1))];
    R2 = [Xv(dTri(:,2)) Yv(dTri(:,2))];
    R3 = [Xv(dTri(:,3)) Yv(dTri(:,3))];

    Rp = [Xp(:) Yp(:)];

    nPoints = numel(Xp);

    % Calculate perpendicular vectors, assuming CCW rotation
    V1 = [-(R2(:,2) - R1(:,2)) (R2(:,1)-R1(:,1))];
    V2 = [-(R3(:,2) - R2(:,2)) (R3(:,1)-R2(:,1))];
    V3 = [-(R1(:,2) - R3(:,2)) (R1(:,1)-R3(:,1))];

    % Calculate number of elements allowed in each iteration
    nPointsInChunk = ceil(maxChunkSizeMB/(numel(V1)*8/(1024^2)));

    % Initialize result vector with NaNs. NaNs will be considered as "out of mesh"
    triIdx = ones([numel(Xp) 1])*NaN;

    % Iterate through points, decide which point is in which triangle
    startIdx = 1;
    msgHdl = waitbar((startIdx-1)/nPoints,'Placing points in mesh');
    while(startIdx < nPoints)
        waitbar((startIdx-1)/nPoints,msgHdl);
        % Tear a chunk out of the point vector
        endIdx = min((startIdx + nPointsInChunk - 1),nPoints);
        cRp = Rp(startIdx:endIdx,:);

        % Check if point is below or above ALL of the lines.
        % (Rp-Ri)*Vi
        d1 = bsxfun(@times,V1(:,1),bsxfun(@plus,cRp(:,1).',-R1(:,1))) + ...
           bsxfun(@times,V1(:,2),bsxfun(@plus,cRp(:,2).',-R1(:,2)));
        d2 = bsxfun(@times,V2(:,1),bsxfun(@plus,cRp(:,1).',-R2(:,1))) + ...
           bsxfun(@times,V2(:,2),bsxfun(@plus,cRp(:,2).',-R2(:,2)));
        d3 = bsxfun(@times,V3(:,1),bsxfun(@plus,cRp(:,1).',-R3(:,1))) + ...
           bsxfun(@times,V3(:,2),bsxfun(@plus,cRp(:,2).',-R3(:,2)));
        isInTri = ((d1 >= 0) & (d2 >= 0) & (d3 >= 0)) | ...
                ((d1 <= 0) & (d2 <= 0) & (d3 <= 0));

        % Find which of these withstood all of the conditions
        [whichTri,whichPoint] = find(isInTri);

        % In case that there is non-unique point assignment, assign to one of them
        [whichPoint,rebuildUnique,~] = unique(whichPoint);
        whichTri = whichTri(rebuildUnique);

        % Store results
        triIdx(whichPoint + startIdx - 1) = whichTri;

        % Advance index
        startIdx = endIdx + 1;

    end

    close(msgHdl);

end


