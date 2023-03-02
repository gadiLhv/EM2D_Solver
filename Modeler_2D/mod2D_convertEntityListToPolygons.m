function [pols,ignoredEnts] = mod2D_convertEntityListToPolygons(dxfEnts)
    % Inputs:
    % 1. dxfEnts - Entity list read by 'importDXF'
    % Outputs:
    % 1. pols - Polygon list
    % 2. ignoredEnts - Entities ignored while converting
    % Initialize output polygon list
    pols = {};
    polCtr = 1;

%    waitHdl = waitbar(0,sprintf('Composed %d polygons from %d/%d entities',numel(pols),0,numel(dxfEnts)));
    ignoredEnts = [];
    % Initialize first polygon starting point
    x = [];
    y = [];
    % Iterate through dxf Entities and start building polygons
    for entIdx = 1:numel(dxfEnts)
        cEntity = dxfEnts{entIdx};

        isLine = strcmp(cEntity.type,'LINE');

        lineIsPoint = false;
        if isLine
            lineIsPoint = ~~min(cEntity.startPoint == cEntity.endPoint,[],2);
        end

        % Currently handles only 'LINE'
        if ~isLine || lineIsPoint
            ignoredEnts = [ignoredEnts ; ...
                    struct('type',cEntity.type,...
                           'number',entIdx)];
            continue;
        end


        if isempty(x)
            % If this is the first point, place both points in
            % array.
            x = [cEntity.startPoint(1) ; cEntity.endPoint(1)];
            y = [cEntity.startPoint(2) ; cEntity.endPoint(2)];
        else
            % If this isn't the first iteration, make a small
            % verification that the previous endpoint is the
            % current start.
            prevEndIsCurrStart = (x(end) == cEntity.startPoint(1)) & ...
                                 (y(end) == cEntity.startPoint(2));

            if prevEndIsCurrStart
                % If indeed this is a concatenated line,
                % store the next endpoint
                x = [x ; cEntity.endPoint(1)];
                y = [y ; cEntity.endPoint(2)];
            else
                % However, if this is not, clear this polygon
                % and throw a warning for a non-closed polygon
                pols(polCtr) = mod2D_createPolygonStruct(x,y);
                polCtr = polCtr + 1;
                fprintf(1,'Polygon #%d is not closed\n',polCtr);
                % And store the current point
                x = [cEntity.startPoint(1) ; cEntity.endPoint(1)];
                y = [cEntity.startPoint(2) ; cEntity.endPoint(2)];
                continue;
            end

            isClosedPolygon = (x(1) == x(end)) & (y(1) == y(end));
            if isClosedPolygon
                x = x(1:(end-1));
                y = y(1:(end-1));

                % Weird bug of closed polygons that have only two points
                if numel(x) == 2
                  fprintf(1,'Polygon dumped due to being a line\n');
                  x = [];
                  y = [];
                  continue;
                end

                % If the polygon is now closed, store and clear
                pols(polCtr) = mod2D_createPolygonStruct(x,y);
                polCtr = polCtr + 1;
                x = [];
                y = [];
            end
        end
%        waitHdl = waitbar(entIdx/numel(dxfEnts),waitHdl,sprintf('Composed %d polygons from %d/%d entities',numel(pols),entIdx,numel(dxfEnts)));
    end

    % If polygon ended and isn't closed, close it manually
    if ~isempty(x)
       pols(polCtr) = mod2D_createPolygonStruct(x,y);
       if numel(x) == 2
          fprintf(1,'Polygon dumped due to being a line\n');
       end
       % If the polygon is now closed, store and clear
       pols(polCtr) = mod2D_createPolygonStruct(x,y);
       polCtr = polCtr + 1;
    end
end
