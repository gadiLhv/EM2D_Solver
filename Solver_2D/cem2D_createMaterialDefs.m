function materialDef = cem2D_createMaterialDefs(varargin)

    % Set default options
    materialDef = struct(...
      'name','default',...  % Material name
      'er',1,...            % Relative permittivity
      'mr',1,...            % Relative permiability
      'type','normal',...   % Could be 'normal', 'metal', 'PEC' or 'PMC'
      'cond_e',1e20,...     % Electric conductivity - used only in case of metallic. Given in Ohm\m
      'cond_m',1e20,...     % Magnetic conductivity - used only in case of metallic. Given in ??? (TBD)
      'tand_e',0,...        % Dielectric loss tangent - used only in case of normal
      'tand_m',0,...        % Magnetic loss tangent - used only in case of normal
      'f_m',0 ...           % The frequency in which the data is measured
      );

    materialDef = misc_validatePropStruct(materialDef,varargin);

end

hold(axHdl,'on');
patch('faces',meshData.tria(meshData.tnum == 1,1:3),'vertices',meshData.vert, ...
    'facecolor','none', ...
    'edgecolor',[0,0,0]) ;