function pols = mod2D_reducePolygonResolution(pols,minDist)
    
    totalVertsBegin = 0;
    totalVertsEnd = 0;
    % Iterate through polygons
    for polIdx = 1:numel(pols)
        cPol = pols{polIdx};
        
        totalVertsBegin = totalVertsBegin + numel(cPol.x);
        cVert = 1;
        
        % Verify that this is necessary here
        distVect = sqrt((cPol.x(2:end) - cPol.x(1:end-1)).^2 + ...
                        (cPol.y(2:end) - cPol.y(1:end-1)).^2);
                     
        if ~sum(distVect < minDist)
            totalVertsEnd = totalVertsEnd + numel(cPol.x);
            continue;
        end
        while cVert < numel(cPol.x)
            nDist = sqrt((cPol.x(cVert + 1) - cPol.x(cVert)).^2 + ...
                         (cPol.y(cVert + 1) - cPol.y(cVert)).^2);
            
            % If distance is to small, clear this element
            if (nDist < minDist)
                cPol.x(cVert + 1) = [];
                cPol.y(cVert + 1) = [];
            end
            
            cVert = cVert + 1;
        end
        
        % Check that last point wasn't removed
        if ~((cPol.x(1) == cPol.x(end)) && (cPol.y(1) == cPol.y(end)))
            cPol.x = [cPol.x ; cPol.x(1)];
            cPol.y = [cPol.y ; cPol.y(1)];
        end
        
        totalVertsEnd = totalVertsEnd + numel(cPol.x);
        
        pols{polIdx} = cPol;
    end
    
    fprintf(1,'Number of vertices reduced by %.2f%%\n',100*(1 - totalVertsEnd/totalVertsBegin));
end