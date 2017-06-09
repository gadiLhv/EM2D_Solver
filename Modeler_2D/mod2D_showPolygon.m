function polHdl = mod2D_showPolygon(axHdl,polStruct,faceColor,edgeColor)
  hold(axHdl)
  
  % Set right the rotation directions (CCW, fill in, CW, fill out)
  [x,y] = polygon2patch(polStruct.x,polStruct.y);
  
  patch (x,y, ...
         'facecolor', faceColor, ...
         'edgecolor', edgeColor);
  
  % Separate the closed loops
%  nanIdxs = find(isnan(polStruct.x));
  polStruct.x = [ polStruct.x ; NaN];
  polStruct.y = [ polStruct.y ; NaN];
  
  nanIdxs = find(isnan(polStruct.x));
  
  % Initialize cell array
  pols = cell([numel(nanIdxs) 1]);
  polIdx = 1;
  while(~isempty(polStruct.x))
    % Find next Nan (indicator for next closed loop)
    nextNan = find(isnan(polStruct.x));
    nextNan = nextNan(1);
    % Read current loop
    pols{polIdx} = [polStruct.x(1:nextNan-1) polStruct.y(1:nextNan-1)];
    polIdx = polIdx + 1;
    % Decimate this loop
    polStruct.x(1:nextNan) = [];
    polStruct.y(1:nextNan) = [];
  end
  
  % Draw the arrow plots
  for polIdx = 1:numel(pols)
    % Read current polygon
    cPol = pols{polIdx};
    
    % Duplicate the first point so the "quiver" will be closed
    cPol_dup = [cPol ; cPol(1,:)];
    % Create arrow direction vector
    cPol_diff = cPol_dup(2:end,:) - cPol_dup(1:end-1,:);
    % Draw the relevant quiver
    h = quiver(cPol(:,1),cPol(:,2),cPol_diff(:,1),cPol_diff(:,2));
    set(h,'color',[0 0 1]);
  end
  
  hold off;
end