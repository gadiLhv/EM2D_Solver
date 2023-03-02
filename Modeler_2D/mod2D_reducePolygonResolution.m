function pols = mod2D_reducePolygonResolution(pols,minDist)

    totalVertsBegin = 0;
    totalVertsEnd = 0;

    dist2vert = @(p1,p2) sqrt(sum((p2 - p1).^2));

    % Empty polygons
    binEmptyPol = false(size(pols));

    % Iterate through polygons
    for polIdx = 1:numel(pols)
        cPol = pols{polIdx};

        % Find loops inside polygon
        loopStart = [1 ; (find(isnan(cPol.x))+1)];
        loopEnd = [(find(isnan(cPol.x))-1) ; numel(cPol.x)];

        totalVertsBegin = totalVertsBegin + numel(cPol.x) - sum(isnan(cPol.x));

        newX = [];
        newY = [];

        % Iterate through loops
        for loopIdx = 1:numel(loopStart)
            newX = [newX ; NaN];
            newY = [newY ; NaN];
            % Iterate through vertices in loop
            cLoop = [cPol.x(loopStart(loopIdx):loopEnd(loopIdx)) cPol.y(loopStart(loopIdx):loopEnd(loopIdx))];
            vertIdx = 1;

            % Check if this is the last loop, which is the only non-self-
            % connected loop.
            isLastLoop = loopIdx == numel(loopStart);
            if(isLastLoop)
                cLoop = [cLoop ; cLoop(1,:)];
            end
            % Iterate through vertices and clean up if necessary
            while vertIdx < size(cLoop,1)
                cVert = cLoop(vertIdx,:);
                nVert = cLoop(vertIdx + 1,:);
                % Check if vertices are too close
                if dist2vert(cVert,nVert) < minDist
                    % If so, remove the next one (arbitrary)
                    cLoop(vertIdx + 1,:) = [];
                else
                    vertIdx = vertIdx + 1;
                end
            end

            % Short sanity check - Check that end is still start
            if sum(cLoop(1,:) == cLoop(end,:)) ~= 2
                cLoop = [cLoop ; cLoop(1,:)];
            end
            totalVertsEnd = totalVertsEnd + size(cLoop,1);

            if(isLastLoop)
                cLoop = cLoop(1:end-1,:);
            end

            % Store new polygon
            newX = [newX ; cLoop(:,1)];
            newY = [newY ; cLoop(:,2)];

        end

        % Remove first NaN
        newX(1) = [];
        newY(1) = [];

        if isempty(newX)
           pols{polIdx} = [];
           binEmptyPol(polIdx) = true;
           continue;
        end

        % And store in a new polygon
        cPol = mod2D_createPolygonStruct(newX,newY);

        pols{polIdx} = cPol;
    end

    fprintf(1,'Number of vertices reduced by %.2f%%\n',100*(1 - totalVertsEnd/totalVertsBegin));
    if sum(binEmptyPol(:))
        fprintf(1,'%d features were deleted\n',sum(binEmptyPol(:)));
        pols(binEmptyPol) = [];
    end
end
