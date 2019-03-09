function parsedVar = mod2D_dxfGroupCodeParser(code,var)
    if (code >= -5) && (code < 0)
        warning('''APP'' codes not supported yet');
        return;
    elseif (code >= 0) && (code < 5)
        % Named strings
        parsedVar = upper(var);
        return;
    elseif code == 5
        % Hexadecimal string of up to 16 characters
        parsedVar = var;
        return;
    elseif code == 6
        % Linetype name
        parsedVar = var;
        return;
    elseif code == 7
        % Text style name
        parsedVar = var;
        return;
    elseif code == 8
        % Layer name
        parsedVar = var;
        return;
    elseif code == 9
        % Variable name identifier
        parsedVar = var;
        return;
    elseif (code >= 10) && (code < 59)
        % Line vertices, points, thicknesses, angles, etc.
        parsedVar = sscanf(var,'%f');
        return;
    elseif code == 60
        % Entity Visibility
        parsedVar = sscanf(var,'%d');
        return;
    elseif code == 62
        % Color number
        parsedVar = sscanf(var,'%d');
        return;
    elseif code == 62
        % Entities follow flag
        parsedVar = sscanf(var,'%d');
        return;
    elseif (code >= 70) && (code < 79)
        % Entities follow flag
        parsedVar = sscanf(var,'%d');
        return;
    elseif (code >= 90) && (code < 100)
        % Entities follow flag
        parsedVar = sscanf(var,'%d');
        return;
    end
end