function pols = mod2D_subtractOvelapingPolygons(pols)
    % Reduce polygon list by subtracting ovelapping polygons
    % This list

    passCtr = 1;

    currentPolIdx = 1;

    binOverlap = false([numel(pols) numel(pols)]);

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

            % In case, for some reason, this returns an empty polygon
            if isempty(newPol)
               % Immediately test the other way around
               newPol = mod2D_booleanOperation(testedPol,cPol,'subtract');
               if isempty(newPol)
%                  % Don't re-store cPol, and remove the tested one
%                  polOverlaps = [polOverlaps ; cTestIdx];
%                  waitHdl = waitbar(currentPolIdx/(numel(pols) - numel(polOverlaps)),waitHdl,sprintf('Intersection check %d/%d',currentPolIdx,(numel(pols) - numel(polOverlaps))));
%                  continue;
                   % Consider this irrecoverable for now.
                   error('Polygons completely overlap');
               end
            else
               % This isn't a full overlap
               if sum(isnan(newPol.x)) > nLoops
                   % Mark that there is an overlap here
                   binOverlap(cTestIdx,currentPolIdx) = true;
                   waitHdl = waitbar(currentPolIdx/(numel(pols) - numel(polOverlaps)),waitHdl,sprintf('Intersection check %d/%d',currentPolIdx,(numel(pols) - numel(polOverlaps))));
%                   % No need to test the other way
                   continue;
               end
            end

            % Test the other way around
            nLoops = sum(isnan(testedPol.x));
            newPol = mod2D_booleanOperation(testedPol,cPol,'subtract');
            if sum(isnan(newPol.x)) > nLoops
                % Mark this as an overlap
                binOverlap(currentPolIdx,cTestIdx) = true;
                waitHdl = waitbar(currentPolIdx/(numel(pols) - numel(polOverlaps)),waitHdl,sprintf('Intersection check %d/%d',currentPolIdx,(numel(pols) - numel(polOverlaps))));
            end
        end

%        % Return clipped polygon to list
%        pols{currentPolIdx} = cPol;
%
%        % Clear polygons from the list
%        pols(polOverlaps) = [];

        currentPolIdx = currentPolIdx + 1;
    end

    % All those that don't overlap

    clippedPolIdxs = [];
    while sum(binOverlap(:))
      % Look at the polygons that have the most overlaps, first.
      [~,sortOrder] = sort(sum(binOverlap,1),'descend');

      cPolIdx = sortOrder(1);

      % Each row that has an odd number of overlaps, is going to be deleted
      cSum = binOverlap(:,cPolIdx).*sum(binOverlap,2);

      % To be deleted
      polsToClip = find(mod(cSum,2));

      % Delete polygons and remove them from table
      cPol = pols{cPolIdx};
      for pIdx = polsToClip(:).'
         clipPol = pols{pIdx};
         cPol = mod2D_booleanOperation(cPol,clipPol,'subtract');

         % In this column,
         binOverlap(binOverlap(:,pIdx),cPolIdx) = false;

         clippedPolIdxs = [clippedPolIdxs ; pIdx];
      end

      % After all polygons were clipped, remove all relevant rows and columns
      binOverlap(polsToClip,:) = [];
      binOverlap(:,polsToClip) = [];

      % Store back polygon
      pols{cPolIdx} = cPol;
    end

    pols(clippedPolIdxs) = [];
    close(waitHdl);
end



