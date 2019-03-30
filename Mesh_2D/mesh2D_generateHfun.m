function fHdl = mesh2D_generateHfun(polygonList,materialAssignment,materialList,meshProps,simProps)
    % While meshing, the mesher sometimes inquires regarding the 
    c0 = physical_constant('speed of light in vacuum');
    e0 = physical_constant('electric constant');
    m0 = physical_constant('mag. constant');

    % TBD: Currently mesh size is set according to the highest simulated frequency.
    f0 = simProps.fMax;

    % Relative mesh edge maximum with respect to wavelength
    relWLmeshMax = meshProps.relWLmeshMax; 

    % Free space wavelength.
    c0 = units('m/sec',[simProps.lengthUnits '*' simProps.freqUnits],c0);

    % Wavelength with respect to current units
    wl0 = c0/f0;
    
    fHdl = @(pts) hfun(pts,wl0,relWLmeshMax,polygonList,materialAssignment,materialList);

end

function maxMeshSize = hfun(pts,wl0,relWLmeshMax,polygonList,materialAssignment,materialList)
    % Default max mesh size
    maxDefualt = wl0*relWLmeshMax;
    
    % Initialize max mesh list
    maxMeshSize = ones([size(pts,1) 1])*maxDefualt;
    % This function actually returns the mesh size
    for polIdx = 1:numel(polygonList)
        cPol = polygonList{polIdx};
        % Obtain material parameters
        cMaterialProps = cem2D_getMaterialPropsFromName(materialAssignment{polIdx},materialList);
        % Get relative permiability and permittivity
        mr = cMaterialProps.mr;
        er = cMaterialProps.er;
        
        % Test all points if they are in this polygon
        binInPoly = isPointInPolygon(pts,[cPol.x cPol.y]);
        
        % Now change the max mesh size, if necessary
        maxMeshSize(binInPoly) =  min(maxMeshSize(binInPoly),wl0*relWLmeshMax/sqrt(er*mr));
    end
end
