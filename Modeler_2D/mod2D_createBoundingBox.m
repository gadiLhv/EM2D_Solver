function bBox = mod2D_createBoundingBox(polList,bBoxAddSpace)
% Returns a polygon structure with the bounding box.

  % Input checks
  if(~exist('polList'))
    errStruct = mod2D_createErrorMsg(...
                  'No polygon list',...
                  'mod2D_calculateBoundingBox');
    throw(errStruct);
  end
  
  if(isempty('polList'))
    errStruct = mod2D_createErrorMsg(...
                  'No polygons in list',...
                  'mod2D_calculateBoundingBox');
    throw(errStruct);
  end
  
  minX = inf;
  maxX = -inf;
  minY = inf;
  maxY = -inf;
  
  % Iterate through polygons
  for cPolCell = polList(:)
    
    cPol = cPolCell{1};
    % Update extrema
    minX = min(minX,min(cPol.x(:)));
    maxX = max(maxX,max(cPol.x(:)));
    minY = min(minY,min(cPol.y(:)));
    maxY = max(maxY,max(cPol.y(:)));
    
  end
   
  bBox = mod2D_createRectangleStruct(...
          [minX minY]+bBoxAddSpace*[-1 1],...
          [maxX maxY]+bBoxAddSpace*[-1 1]);
end