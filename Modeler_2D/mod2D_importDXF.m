function [header,tables,entities] = mod2D_importDXF(fileName)

    % Read entire file
    fHdl = fopen(fileName);
    dxfRows = textscan(fHdl,'%d%s','delimiter','\n');
    fclose(fHdl);
    valType = dxfRows{1};
    dxfCode = dxfRows{2};
    clear('dxfRows');

    % Verify EOF
    if strcmp(dxfCode{end},'EOF')
        valType = valType(1:end-1);
        dxfCode = dxfCode(1:end-1);
    else
        warning('No ''EOF'' indicator found');
    end

    % Disect to sections
%    sectionIdxs = strmatch('SECTION',dxfCode,'exact');
    matchSection = @(s) strcmp(s,'SECTION');
    sectionIdxs = find(cellfun(matchSection,dxfCode));
    binMap = false(size(valType));
    binMap(sectionIdxs) = true;
    binMap = binMap & (valType == 0);
    sectionIdxs = find(binMap);

    header = [];
    tables = [];
    entities = [];
    % Parse sections one by one
    sectBoundaries = [sectionIdxs(:) [sectionIdxs(2:end) ; numel(dxfCode)]];
    for sectCtr = 1:numel(sectionIdxs)
        cSectIdxs = sectBoundaries(sectCtr,:);

        cSect_valType = valType((cSectIdxs(1) + 1):(cSectIdxs(2) - 1));
        cSect_dxfCode = dxfCode((cSectIdxs(1) + 1):(cSectIdxs(2) - 1));

        % Type of section
        sectionType = cSect_dxfCode{1};

        switch sectionType
            case 'HEADER'
                header = header_handler(cSect_valType(2:end),cSect_dxfCode(2:end));
            case 'TABLES'
                tables = tables_handler(cSect_valType(2:end),cSect_dxfCode(2:end));
            case 'ENTITIES'
                entities = entity_handler(cSect_valType(2:end),cSect_dxfCode(2:end));
            case 'BLOCKS'
                fprintf(1,'No handler for ''BLOCKS'' just now.\n');

        end
    end

end

function entities = entity_handler(groupCode,variable)

    % Disect to separate entities
    entityIdxs = find(groupCode == 0);

    entities = cell([numel(entityIdxs)-1 1]);

    for entIdx = 1:(numel(entityIdxs) - 1)

        cCode = groupCode(entityIdxs(entIdx):(entityIdxs(entIdx + 1) - 1));
        cVar = variable(entityIdxs(entIdx):(entityIdxs(entIdx + 1) - 1));

        % First store original type
        cEnt = struct('type',cVar{1});
        % Give graphical interface allocation
        switch cEnt.type
            case 'LINE'
                cEnt.mod2D_assign = 'Graphical';    % TODO: Define other types for text and other stuff?
                cEnt.startPoint = [0 0 0];          % Start point, [X Y Z]
                cEnt.endPoint = [0 0 0];            % End point, [X Y Z]
                cEnt.lineType = [];                 % According to LTYPE tables
                cEnt.layer = 'Default';             % Arbitrary value. Layer assignment is imperative to determine continuity of entities
                cEnt.thickness = 0;
                % Extract all properties
                for coorIdx = 1:3
                    cEnt.startPoint(coorIdx) = searchCodeReturnVal(cCode,cVar,coorIdx*10,cEnt.startPoint(coorIdx));
                    cEnt.endPoint(coorIdx) = searchCodeReturnVal(cCode,cVar,coorIdx*10 + 1,cEnt.endPoint(coorIdx));
                end
                % Misc.
                cEnt.layer = searchCodeReturnVal(cCode,cVar,8,cEnt.layer);
                cEnt.lineType = searchCodeReturnVal(cCode,cVar,6,cEnt.lineType);
                cEnt.thickness = searchCodeReturnVal(cCode,cVar,39,cEnt.thickness);
            case 'ARC'
                cEnt.mod2D_assign = 'Graphical';    % TODO: Define other types for text and other stuff?
                cEnt.startPoint = [1 0 0];          % Start point, [X Y Z]
                cEnt.endPoint = [0 1 0];            % End point, [X Y Z]
                cEnt.centerPoint = [0 0 0];         % Center of arc
                cEnt.startAngle = 0;                % Start angle of rotation
                cEnt.endAngle = 90;                 % End angle. Converted here to radians
                cEnt.radius = 1;                    % Radius of arc
                cEnt.lineType = [];                 % According to table
                cEnt.layer = 'Default';             % Arbitrary value. Layer assignment is imperative to determine continuity of entities
                cEnt.thickness = 0;
                for coorIdx = 1:3
                    cEnt.centerPoint(coorIdx) = searchCodeReturnVal(cCode,cVar,coorIdx*10,cEnt.centerPoint(coorIdx));
                end
                cEnt.radius = searchCodeReturnVal(cCode,cVar,40,cEnt.radius);
                cEnt.startAngle = searchCodeReturnVal(cCode,cVar,50,cEnt.startAngle);
                cEnt.endAngle = searchCodeReturnVal(cCode,cVar,51,cEnt.endAngle);
                % Convert center & start\end angles to actual start-end
                % coordinates (useful for geometrical processing, later).
                cEnt.startPoint = cEnt.centerPoint + cEnt.radius*[cos(cEnt.startAngle*pi/180) sin(cEnt.startAngle*pi/180) 0];
                cEnt.endPoint = cEnt.centerPoint + cEnt.radius*[cos(cEnt.endAngle*pi/180) sin(cEnt.endAngle*pi/180) 0];
                % Misc.
                cEnt.layer = searchCodeReturnVal(cCode,cVar,8,cEnt.layer);
                cEnt.lineType = searchCodeReturnVal(cCode,cVar,6,cEnt.lineType);
                cEnt.thickness = searchCodeReturnVal(cCode,cVar,39,cEnt.thickness);
            case 'CIRCLE'
                cEnt.mod2D_assign = 'Graphical';    % TODO: Define other types for text and other stuff?
                cEnt.startPoint = [0 0 0];          % Start point, [X Y Z]
                cEnt.centerPoint = [0 0 0];         % Center of arc
                cEnt.radius = [];                   % Radius of arc
                cEnt.lineType = [];                 % According to
                cEnt.layer = 'Default';             % Arbitrary value. Layer assignment is imperative to determine continuity of entities
                cEnt.thickness = 0;
                % Read center [X Y Z] coordinates
                for coorIdx = 1:3
                    cEnt.centerPoint(coorIdx) = searchCodeReturnVal(cCode,cVar,coorIdx*10,cEnt.centerPoint(coorIdx));
                end
                % Other circle properties
                cEnt.radius = searchCodeReturnVal(cCode,cVar,40,cEnt.radius);
                % Misc.
                cEnt.layer = searchCodeReturnVal(cCode,cVar,8,cEnt.layer);
                cEnt.lineType = searchCodeReturnVal(cCode,cVar,6,cEnt.lineType);
                cEnt.thickness = searchCodeReturnVal(cCode,cVar,39,cEnt.thickness);
            otherwise
                warning(sprintf('Cannot handle entities of type ''%s'' yet',cVar{1}));
        end
    % Store current entity
    entities{entIdx} = cEnt;

    end

end

function retVal = searchCodeReturnVal(code,var,codeToFind,defaultVal)
    retVal = defaultVal;
    retValIdx = find(code == codeToFind);
    if ~isempty(retValIdx)
        retVal = mod2D_dxfGroupCodeParser(code(retValIdx),var{retValIdx});
    end
end

function tables = tables_handler(groupCode,variable)
    varIdx = 1;
    cTable = [];
    aggTable = 0;
    % The field that must exist is 'Layer'

    isTable = @(cVal) strcmp('TABLE',cVal);
    isTableEnd = @(cVal) strcmp('ENDTAB',cVal);
    % Find all the indexes that define the beginning of a table
    tableStartIdxs = find((groupCode == 0) & cellfun(isTable,variable));
    tableEndIdxs = find((groupCode == 0) & cellfun(isTableEnd,variable));

    ignoreCtr = 0;

    % Initialize the table container
    tables = cell([numel(tableStartIdxs) 1]);
    % Iterate through them
    for tableIdx = 1:numel(tableStartIdxs)
        cCode = groupCode((tableStartIdxs(tableIdx) + 1):(tableEndIdxs(tableIdx) - 1));
        cVar = variable((tableStartIdxs(tableIdx) + 1):(tableEndIdxs(tableIdx) - 1));

        % First group code is the table type.
        % Handle table type
        switch cVar{1}
            case 'LTYPE'
                cCode = cCode(2:end);
                cVar = cVar(2:end);
                % Count how many line types will be defined here
                isLtype = @(cVal) strcmp('LTYPE',cVal);
                lTypesIdxs = find((cCode == 0) & cellfun(isLtype,cVar));

                cTables = [];

                lTypesIdxs = [lTypesIdxs ; numel(cCode)+1];
                % Iterate through ltypes and fill structures
                for lIdx = 1:(numel(lTypesIdxs) - 1)
                    % Snip out relevant parts
                    LtypeCode = cCode((lTypesIdxs(lIdx)+1):lTypesIdxs(lIdx+1)-1);
                    LtypeVar = cVar((lTypesIdxs(lIdx)+1):lTypesIdxs(lIdx+1)-1);

                    % Structre template
                    cTable = struct('tableType','LTYPE',...
                                    'name',[],...
                                    'description',[],...
                                    'handle',[],...
                                    'alignment',[],...
                                    'numberOfElements',[],...
                                    'totalPatternLength',[],...
                                    'shapeLength',[]);

                    % Iterate through fields and fill one by one
                    for fieldIdx = 1:numel(LtypeCode)
                        cLtypeCode = LtypeCode(fieldIdx);
                        cLtypeVar = mod2D_dxfGroupCodeParser(cLtypeCode,LtypeVar{fieldIdx});
                        switch cLtypeCode
                            case 2
                                cTable.name = cLtypeVar;
                            case 3
                                cTable.description = cLtypeVar;
                            case 5
                                cTable.handle = cLtypeVar;
                            case 40
                                cTable.totalPatternLength = cLtypeVar;
                            case 49
                                cTable.shapeLength = [cTable.shapeLength ; cLtypeVar];
                            case 72
                                cTable.alignment = cLtypeVar;
                            case 73
                                cTable.numberOfElements = cLtypeVar;
                            otherwise
                                ignoreCtr = ignoreCtr + 1;

                        end
                    end
                    % Store current table
                    cTables = [cTables ; cTable];
                end

                % Store in return value, the table cell array
                tables{tableIdx} = cTables;
            case 'LAYER'
                cCode = cCode(2:end);
                cVar = cVar(2:end);
                % Count how many line types will be defined here
                isLayer = @(cVal) strcmp('LAYER',cVal);
                layerIdxs = find((cCode == 0) & cellfun(isLayer,cVar));

                cLayers = [];

                layerIdxs = [layerIdxs ; numel(cCode)+1];
                % Iterate through layers and fill structures
                for lIdx = 1:(numel(layerIdxs) - 1)
                    % Snip out relevant parts
                    layerCode = cCode((layerIdxs(lIdx)+1):layerIdxs(lIdx+1)-1);
                    layerVar = cVar((layerIdxs(lIdx)+1):layerIdxs(lIdx+1)-1);

                    % Structre template
                    cLayer = struct('name',[],...
                                    'colorNumber',[],...
                                    'lineTypeName',[],...
                                    'plottingFlag',[],...
                                    'lineWeight',[],...
                                    'hardPtr',[]);

                    for fieldIdx = 1:numel(layerCode)
                        cLayerCode = layerCode(fieldIdx);
                        cLayerVar = mod2D_dxfGroupCodeParser(cLayerCode,layerVar{fieldIdx});

                        switch cLayerCode
                            case 2
                                cLayer.name = cLayerVar;
                            case 6
                                cLayer.lineTypeName = cLayerVar;
                            case 62
                                cLayer.colorNumber = cLayerVar;
                            case 290
                                cLayer.plottingFlag = cLayerVar;
                            case 370
                                cLayer.lineWeight = cLayerVar;
                            case 390
                                cLayer.lineWeight = cLayerVar;
                            otherwise
                                ignoreCtr = ignoreCtr + 1;
                        end
                    end
                    % Store current table
                    cLayers = [cLayers ; cLayer];
                end

                fprintf(1,'%d ''LAYER'' codes ignored\n');
                % Store in output list
                tables{tableIdx} = cLayers;
            case 'STYLE'
                cCode = cCode(2:end);
                cVar = cVar(2:end);
                % Count how many line types will be defined here
                isStyle = @(cVal) strcmp('STYLE',cVal);
                styleIdxs = find((cCode == 0) & cellfun(isStyle,cVar));

                cStyles = [];

                styleIdxs = [styleIdxs ; numel(cCode)+1];
                % Iterate through ltypes and fill structures
                for sIdx = 1:(numel(styleIdxs) - 1)
                    % Snip out relevant parts
                    styleCode = cCode((styleIdxs(sIdx)+1):styleIdxs(sIdx+1)-1);
                    styleVar = cVar((styleIdxs(sIdx)+1):styleIdxs(sIdx+1)-1);

                    % Structre template
                    cStyle = struct('name',[],...
                                    'primaryFontName',[],...
                                    'bigFontName',[],...
                                    'fixedTextHeight',[],...
                                    'widthFactor',[],...
                                    'lastHeightUsed',[],...
                                    'obliqueAngle',[],...
                                    'genFlag',[]);

                    for fieldIdx = 1:numel(styleCode)
                        cStyleCode = styleCode(fieldIdx);
                        cStyleVar = mod2D_dxfGroupCodeParser(cStyleCode,styleVar{fieldIdx});

                        switch cStyleCode
                            case 2
                                cStyle.name = cStyleVar;
                            case 3
                                cStyle.primaryFontName = cStyleVar;
                            case 4
                                cStyle.bigFontName = cStyleVar;
                            case 40
                                cStyle.fixedTextHeight = cStyleVar;
                            case 41
                                cStyle.widthFactor = cStyleVar;
                            case 42
                                cStyle.lastHeightUsed = cStyleVar;
                            case 50
                                cStyle.obliqueAngle = cStyleVar;
                            case 71
                                cStyle.genFlag = cStyleVar;
                            otherwise
                                ignoreCtr = ignoreCtr + 1;
                        end
                    end
                    % Add to style list
                    cStyles = [cStyles ; cStyle];
                end

                % Store in output list
                tables{tableIdx} = cStyles;
            otherwise
                warning(sprintf('Table type ''%s'' is currently not supported',cVar{1}));
        end
    end


    if ignoreCtr
        fprintf(1,'Total of %d codes ignored while reading DXF\n',ignoreCtr);
    end
end

function header = header_handler(groupCode,variable)
    % Iterate until list is done
    varIdx = 1;
    while varIdx <= numel(variable)
        cCode = groupCode(varIdx);
        cVar = variable{varIdx};

        % Check if file ended
        if strcmp(cVar,'ENDSEC')
            return;
        end

        % Parse out the variable value in the correct format
        cVal = mod2D_dxfGroupCodeParser(groupCode(varIdx + 1),variable{varIdx + 1});
        eval(['header.' cVar(2:end) ' = cVal;']);

        varIdx = varIdx + 2;
    end


end

