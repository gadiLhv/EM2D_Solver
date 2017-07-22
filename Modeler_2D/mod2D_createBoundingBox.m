function bBox = mod2D_createBoundingBox(polList,meshProps,simProps,bgMaterial)
% Returns a polygon structure with the bounding box.

c0 = physical_constant('speed of light in vacuum');
e0 = physical_constant('electric constant');
m0 = physical_constant('mag. constant');

% Calculate wavelength in current units:
% 1. Choose maximum frequency as reference for shortest simulated wavelength
fMin = simProps.fMin;
% 2. Convert light speed to '(length units)*(frequency units)'
c0 = units('m/sec',[simProps.lengthUnits '*' simProps.freqUnits],c0);
% 3. Wavelength in current units is given by frequency with current units
wl = c0/fMin;

er = bgMaterial.er;
mr = bgMaterial.mr;
bBoxAddSpace = wl*meshProps.boundingBoxAddSpace/sqrt(er*mr);

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
for cPolCell = polList(:).'
  
  cPol = cPolCell{1};
  % Update extrema
  minX = min(minX,min(cPol.x(:)));
  maxX = max(maxX,max(cPol.x(:)));
  minY = min(minY,min(cPol.y(:)));
  maxY = max(maxY,max(cPol.y(:)));
  
end
 
bBox = mod2D_createRectangleStruct(...
        [minX minY]+bBoxAddSpace*[-1 -1],...
        [maxX maxY]+bBoxAddSpace*[1 1]);
        
end