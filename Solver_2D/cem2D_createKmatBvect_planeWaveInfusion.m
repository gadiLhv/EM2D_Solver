function cem2D_createKmatBvect_planeWaveInfusion(meshData,materialList,materialAssign,simProps,f_sim,pwStruct)
% Inputs: 
% 1. meshData - Structure with all of the mesh data
% 2. materialList - Cell array with all the material data structures.
% 3. materialAssign - Material assignements per face (object).
% 4. simProps - General simulation properties
% 5. f_sim - Frequency (in simulation units) of current solution

distTH = 1e-10;

c0 = physical_constant('speed of light in vacuum');
e0 = physical_constant('electric constant');
m0 = physical_constant('mag. constant');

% Convert light speed to correct units
c0 = units('m/sec',[simProps.lengthUnits '*' simProps.freqUnits],c0);
% Calculate wave number
k0 = 2*pi*f_sim/c0;

% Extract background material (allways first face) properties
cMaterialProps = cem2D_getMaterialPropsFromName(materialAssign{1},materialList);

nVerts = size(meshData.vert,1);

% Initial FEM matrix and source vector
K = zeros([1 1]*nVerts);
b = zeros([nVerts 1]);

% Get relative permiability and permittivity
mr = cMaterialProps.mr;
er = cMaterialProps.er;

% Recover all triangles relevant to background face
triBinMap = meshData.tnum == 1;
% Store all triplets 
triTriplets = meshData.tria(triBinMap,:);

% Find edge node pairs in this set of triangles.
edgePairs = cem2D_findEdgeNodes(triTriplets,meshData.etri);

% Extract only edges on the outermost part
% 1. Find the extrema of the bounding box
xRange = [min(meshData.vert(:,1)) max(meshData.vert(:,1))];
yRange = [min(meshData.vert(:,2)) max(meshData.vert(:,2))];
% 2. Extract all pairs
Rpairs1 = meshData.vert(edgePairs(:,1),:);
Rpairs2 = meshData.vert(edgePairs(:,2),:);
% 3. Determine if the pairs are on any of the extemum. Decide with threshold, 
%    rather than simply comparing. Mesh isn't perfect...
isOnEdge =  (inTH(Rpairs1(:,1),xRאיילוןange(1),distTH) & inTH(Rpairs2(:,1),xRange(1),distTH)) | ...
            (inTH(Rpairs1(:,1),xRange(2),distTH) & inTH(Rpairs2(:,1),xRange(2),distTH)) | ...
            (inTH(Rpairs1(:,2),yRange(1),distTH) & inTH(Rpairs2(:,2),yRange(1),distTH)) | ...
            (inTH(Rpairs1(:,2),yRange(2),distTH) & inTH(Rpairs2(:,2),yRange(2),distTH));
% 4. Extract only outer edge pairscem2D_findEdgeNodes
edgePairs = edgePairs(isOnEdge,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debug: paint triangles and faces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  figure;
%  patch('faces',triTriplets,...
%        'vertices',meshData.vert, ...
%        'facecolor',[1.,1.,1.], ...
%        'edgecolor',[0,0,0]) ;
%  hold on; 
%  axis image off;
%  edgesX = [meshData.vert(edgePairs(:,1),1) meshData.vert(edgePairs(:,2),1)].';
%  edgesY = [meshData.vert(edgePairs(:,1),2) meshData.vert(edgePairs(:,2),2)].';
%  plot(edgesX,edgesY,'-','linewidth',3);
%  hold off;       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debug: paint triangles and faces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate absoulte distance from origin
rho1 = meshData.vert(edgePairs(:,1),:);
rho2 = meshData.vert(edgePairs(:,2),:);

% Determine tangent vector (phi)
dr = rho2 - rho1);
drAbs = sqrt(sum(dr.^2,2));
phi = dr./drAbs;

kappa = 1./sqrt(sum(rho1.^2,2));

% Determine normal (outer) by performing cross between the Z direction and the 
% tangent vector.
n = [-phi(:,2) phi(:,1)];

% In case of plane wave infusion, it is only necessary to calculate the "q" 
% coefficients, to insert into the 'b' vector.

% Calculate segment lengths
l_s = sqrt(sum((rho2 - rho1).^2,2));

% Calculate explicit propagation direction
kx = cos(pwStruct.phiPropagation*pi/180);
ky = sin(pwStruct.phiPropagation*pi/180);


% Freespace impedance
eta = sqrt(m0*mr./(e0*er));

% Derivate the incident field in both directions
dphidn = -1i*k0*(kx*n(:,1) + ky*n(:,2)).*exp(-1i*k0*(kx*rho(:,1) + ky*rho(:,2)));
d2phidphi2 = -(k0^2)*((kx*n(:,1) + ky*n(:,2))^2).*exp(-1i*k0*(kx*rho(:,1) + ky*rho(:,2)));
% Build q vector
q = dphidn+ ...
    [1i*k0 + 0.5*kappa - 1i*kappa.^2./(8*(1i*kappa - k0))) - ...
    (0.5*1i./(kappa - k0)).*d2phidphi2;

    
switch simProps.polarizationType
  case 'TE'
    % First 
    q = (pwStruct.amp./eta)*q;
  case 'TM'
    q = (pwStruct.amp)*q;
  case 'TEM'
    error('There is no such thing as a plane wave in ''TEM''');
end
  
% No need to add elements to K matrix because this imposes only sources. 
% The radiation boundary condtions are handled separately. 

% To the sources vector, however, here it is!
b(edgePairs(:,1)) = b(edgePairs(:,1)) + q;


end

function binAns = inTH(a,b,TH)
binAns = abs(a - b) < TH;
end