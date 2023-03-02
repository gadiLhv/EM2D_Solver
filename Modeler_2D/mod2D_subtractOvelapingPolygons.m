function pols = mod2D_subtractOvelapingPolygons(pols)
    % Reduce polygon list by subtracting ovelapping polygons
    % This list

    passCtr = 1;

    currentPolIdx = 1;

    waitHdl = waitbar(0,sprintf('Intersection check %d/%d',currentPolIdx,numel(pols)));
    while currentPolIdx < numel(pols)

        % Check 1st item in the list
        cPol = pols{currentPolIdx};

        % Count overlaps with this polygon
        polOverlaps = [];

        waitHdl = waitbar(currentPolIdx/(numel(pols) - numel(polOverlaps)),waitHdl,sprintf('Intersection check %d/%d',currentPolIdx,(numel(pols) - numel(polOverlaps))));

        % Check overlap of current polygon with polygons
        % below it in the list
        for cTestIdx = (currentPolIdx + 1):numel(pols)

            % Check both directions of overlap
            testedPol = pols{cTestIdx};

            % Test both ways
            % Check previous number of loops
            nLoops = sum(isnan(cPol.x));
            % Attempt to clip the polygons one way
            newPol = mod2D_booleanOperation(cPol,testedPol,'subtract');
            if sum(isnan(newPol.x)) > nLoops
                cPol = newPol;
                % Count this overlap, and consider polygons to remove
                polOverlaps = [polOverlaps ; cTestIdx];
                waitHdl = waitbar(currentPolIdx/(numel(pols) - numel(polOverlaps)),waitHdl,sprintf('Intersection check %d/%d',currentPolIdx,(numel(pols) - numel(polOverlaps))));
                % No need to test the other way
                continue;
            end

            % Test the other way around
            nLoops = sum(isnan(testedPol.x));
            newPol = mod2D_booleanOperation(testedPol,cPol,'subtract');
            if sum(isnan(newPol.x)) > nLoops
                cPol = newPol;
                % If this is the case, the new polygon is re-instated as the original one.
                polOverlaps = [polOverlaps ; cTestIdx];
                waitHdl = waitbar(currentPolIdx/(numel(pols) - numel(polOverlaps)),waitHdl,sprintf('Intersection check %d/%d',currentPolIdx,(numel(pols) - numel(polOverlaps))));
            end
        end

        % Return clipped polygon to list
        pols{currentPolIdx} = cPol;

        % Clear polygons from the list
        pols(polOverlaps) = [];

        currentPolIdx = currentPolIdx + 1;
    end

%    % Split polygons if necessary
%    polCtr = 1;
%    while polCtr <= numel(pols)
%       cPol = pols{polCtr};
%
%    end

    close(waitHdl);
end
