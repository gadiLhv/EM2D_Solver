function [A,B] = cem2D_createABmat_lossless(meshData,meshProps,materialList,materialAssign,simProps,f_sim)
% Create a matrix pair for eigenvalue decomposition - Port mode analysis
%
% Inputs:
% 1. meshData - Structure with all of the mesh data
% 2. materialList - Cell array with all the material data structures.
% 3. materialAssign - Material assignements per face (object).
% 4. simProps - General simulation properties
% 5. f_sim - Frequency (in simulation units) of current solution

c0 = physical_constant('speed of light in vacuum');
e0 = physical_constant('electric constant');
m0 = physical_constant('mag. constant');

% Convert light speed to correct units
f_sim = units(simProps.freqUnits,'Hz',f_sim);
% Calculate wave number
k0 = 2*pi*f_sim/c0;

nEdges = size(meshData.edge,1);

% Initial FEM matrix and source vector
A = zeros([1 1]*nEdges);
B = A;

% Build matrix face by face
for faceIdx = 1:numel(meshData.face)
    % Extract current parts material properties
    cMaterialProps = cem2D_getMaterialPropsFromName(materialAssign{faceIdx},materialList);

    % If this specific face is metal, continue
    if(~strcmp(cMaterialProps.type,'normal'))
    continue;

    % Otherwise, start building the cell

    % Get relative permiability and permittivity
    mr = cMaterialProps.mr;
    er = cMaterialProps.er;

    imr = 1/mr;

    % Recover all triangles relevant to current face
    triBinMap = meshData.tnum == faceIdx;
    % Store all triplets
    triTriplets = meshData.tria(triBinMap,:);

    % Calculate edge lengths
    eL = cem_calcEdgeLengths(vert_m,triTriplets);

    % Create interpolants for elements
    vert_m = units(simProps.lengthUnits,'m',meshData.vert);
    [a,b,c,Det] = cem2D_createInterpolantCoeffs1st(vert_m,triTriplets);

    % Initialize triangle element sub-matrices
    Ae = imr*Ee(eL,Det) + (er*k0^2)*Fe(eL,Det,b,c);

    Be_tt = imr*Fe(eL,Det,b,c);
    Be_tz =



end % For main function, "cem2D_createABmat_lossless"

%
%
%    % Initialize triangle element sub-matrices
%    Ke = zeros([3 3 size(c,1)]);
%    for i = 1:3
%        for j = 1:3
%            Ke(i,j,:) = ((4*Det).^(-1)).*(ax.*b(:,i).*b(:,j) + ...
%            ay.*c(:,i).*c(:,j)) + (Det/12).*beta.*(1 + (i == j));
%        end
%    end
%
%    % Add the elements from the local triangle element sub matrices
%    for e = 1:size(Ke,3)
%        vIdxs = [triTriplets(e,1) triTriplets(e,2) triTriplets(e,3)];
%        K(vIdxs,vIdxs) = K(vIdxs,vIdxs) + Ke(:,:,e);
%    end
%
%end
%
%b = zeros([size(K,1) 1]);

end

function cE = Ee(edgeL,detE)

edgeL = permute(edgeL,[2 3 1]);
detE = permute(detE,[3 2 1]);

% Calculate the Eij matrix per element
cE = bsxfun(@times,bsxfun(@times,edgeL,permute(edgeL,[2 1 3])),1./detE);

end

function cF = Fe(edgeL,detE,b,c)

edgeL = permute(edgeL,[2 3 1]);
detE = permute(detE,[3 2 1]);

b = permute(b,[2 3 1]);
c = permute(c,[2 3 1]);

f = @(i,j) b(i,:,:).*b(j,:,:) + c(i,:,:).*c(j,:,:);

Fe11 = (edgeL(1,:,:).*edgeL(1,:,:)./(24*detE)).*(f(2,2) - f(1,2) + f(1,1));
Fe12 = (edgeL(1,:,:).*edgeL(2,:,:)./(48*detE)).*(f(2,3) - f(2,2) - 2*f(1,3) + f(1,2));
Fe13 = (edgeL(1,:,:).*edgeL(3,:,:)./(48*detE)).*(f(2,1) - 2*f(2,3) - f(1,1) + f(1,3));

Fe22 = (edgeL(2,:,:).*edgeL(2,:,:)./(24*detE)).*(f(3,3) - f(2,3) + f(2,2));
Fe23 = (edgeL(2,:,:).*edgeL(3,:,:)./(48*detE)).*(f(3,1) - f(3,3) - 2*f(2,1) + f(2,3));

Fe33 = (edgeL(3,:,:).*edgeL(3,:,:)./(24*detE)).*(f(1,1) - f(1,3) + f(3,3));

cF = [Fe11 Fe12 Fe13 ; Fe12 Fe22 Fe23 ; Fe13 Fe23 Fe33];

end

