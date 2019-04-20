function polStruct = mod2D_createPolygonStruct(X,Y)
    % Returns a structure with the (x,y) coordinates of the 
    % Fields:
    % x - X coordinates, aligned in a single column
    % y - Y coordinates, aligned in a single column
    % nParts - Number of closed loops that composes polygons
    % structType - string. Name of the element type
    %
    % No need to repeat first coordinate twice

    X = X(:);
    Y = Y(:);
    
    % Count loops
    nanIdxs = find(isnan(X));
    nLoops = numel(nanIdxs) + 1;
    
    loopStart = [1 (nanIdxs(:) + 1).'];
    loopEnd = [(nanIdxs(:) - 1).' numel(X)];
    
    % Verify that last loop is not connected to itself
    loopStart = loopStart(end);
    loopEnd = loopEnd(end);
    
    if (X(loopStart) == X(loopEnd)) && (Y(loopStart) == Y(loopEnd))
        warning('Last loop of polygon cannot be self connected. Removing...');
        X = X(1:(end-1));
        Y = Y(1:(end-1));
    end
    polStruct = struct('x',X(:),'y',Y(:),'nParts',sum(isnan(X(:))) + 1,'structType','Polygon');

end