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

    vert_m = units(simProps.lengthUnits,'m',meshData.vert);

    [eIdxs,nEdges,~] = mesh2D_createEdgeIndexing(meshData);

    % Initial FEM matrix and source vector
    A = zeros([1 1]*nEdges);
    B_tt = A;
    B_tz = A;
    B_zt = A;
    B_zz = A;

    % Build matrix face by face
    for faceIdx = 1:numel(meshData.face)
        % Extract current parts material properties
        cMaterialProps = cem2D_getMaterialPropsFromName(materialAssign{faceIdx},materialList);

        % If this specific face is metal, continue
        if(~strcmp(cMaterialProps.type,'normal'))
            continue;
        end

        % Otherwise, start building the cell

        % Get relative permiability and permittivity
        mr = cMaterialProps.mr;
        er = cMaterialProps.er;

        imr = 1/mr;

        % Recover all triangles relevant to current face
        triBinMap = meshData.tnum == faceIdx;
        % Store all triplets
        triTriplets = meshData.tria(triBinMap,:);
        edgeTripliets = eIdxs(triBinMap,:);

        % Calculate edge lengths
        eL = cem_calcEdgeLengths(vert_m,triTriplets);

        % Create interpolants for elements
        [a,b,c,Det] = cem2D_createInterpolantCoeffs1st(vert_m,triTriplets);

        % A matrix is simple enough
        Ae = imr*Ee(eL,Det) + (er*k0^2)*Fe(eL,Det,b,c);

        % B matrix top left sub-matrix is still pretty simple
        Be_tt = imr*Fe(eL,Det,b,c);

        % Here it gets nasty
        eL = permute(eL,[2 3 1]);
        Det = permute(Det,[3 2 1]);
        b = permute(b,[2 3 1]);
        c = permute(c,[2 3 1]);

        % "f" elements defined in page 242, "Finite element method in electromagnetics", 1st edition.
        f = @(i,j) b(i,:,:).*b(j,:,:) + c(i,:,:).*c(j,:,:);

        % Trans-pol matrices
        Be_tz_11 = (imr*eL(1,:,:)./(Det*12)).*(f(2,1) - f(1,1));
        Be_tz_12 = (imr*eL(1,:,:)./(Det*12)).*(f(2,2) - f(1,2));
        Be_tz_13 = (imr*eL(1,:,:)./(Det*12)).*(f(2,3) - f(1,3));

        Be_tz_21 = (imr*eL(2,:,:)./(Det*12)).*(f(3,1) - f(2,1));
        Be_tz_22 = (imr*eL(2,:,:)./(Det*12)).*(f(3,2) - f(2,2));
        Be_tz_23 = (imr*eL(2,:,:)./(Det*12)).*(f(3,3) - f(2,3));

        Be_tz_31 = (imr*eL(3,:,:)./(Det*12)).*(f(1,1) - f(3,1));
        Be_tz_32 = (imr*eL(3,:,:)./(Det*12)).*(f(1,2) - f(3,2));
        Be_tz_33 = (imr*eL(3,:,:)./(Det*12)).*(f(1,3) - f(3,3));

        Be_tz = [   Be_tz_11 Be_tz_12 Be_tz_13 ; ...
                    Be_tz_21 Be_tz_22 Be_tz_23 ; ...
                    Be_tz_31 Be_tz_32 Be_tz_33];

        Be_zt_11 = (imr*eL(1,:,:)./(Det*12)).*(f(1,2) - f(1,1));
        Be_zt_12 = (imr*eL(1,:,:)./(Det*12)).*(f(1,3) - f(1,2));
        Be_zt_13 = (imr*eL(1,:,:)./(Det*12)).*(f(1,1) - f(1,3));

        Be_zt_21 = (imr*eL(2,:,:)./(Det*12)).*(f(2,2) - f(2,1));
        Be_zt_22 = (imr*eL(2,:,:)./(Det*12)).*(f(2,3) - f(2,2));
        Be_zt_23 = (imr*eL(2,:,:)./(Det*12)).*(f(2,1) - f(2,3));

        Be_zt_31 = (imr*eL(3,:,:)./(Det*12)).*(f(3,2) - f(3,1));
        Be_zt_32 = (imr*eL(3,:,:)./(Det*12)).*(f(3,3) - f(3,2));
        Be_zt_33 = (imr*eL(3,:,:)./(Det*12)).*(f(3,1) - f(3,3));

        Be_zt = [   Be_zt_11 Be_zt_12 Be_zt_13 ; ...
                    Be_zt_21 Be_zt_22 Be_zt_23 ; ...
                    Be_zt_31 Be_tz_32 Be_zt_33];

        % And finally, cross-pol elements:
        Be_zz = zeros([3 3 numel(Det)]);
        for i = 1:3
            for j = 1:3
                Be_zz(i,j,:) = f(i,j)./(4*Det) - (k0^2)*er*(1 + (i == j)).*Det/12;
            end
        end

        % Iteratively add sub-matrices
        for e = 1:size(Ae,3)
            ceIdxs = edgeTripliets(e,:);
            A(ceIdxs,ceIdxs) = A(ceIdxs,ceIdxs) + Ae(:,:,e);
            B_tt(ceIdxs,ceIdxs) = B_tt(ceIdxs,ceIdxs) + Be_tt(:,:,e);
            B_tz(ceIdxs,ceIdxs) = B_tz(ceIdxs,ceIdxs) + Be_tz(:,:,e);
            B_zt(ceIdxs,ceIdxs) = B_zt(ceIdxs,ceIdxs) + Be_zt(:,:,e);
            B_zz(ceIdxs,ceIdxs) = B_zz(ceIdxs,ceIdxs) + Be_zz(:,:,e);
        end

    end % Per-face for loop

    Zz = zeros(size(A));

%    A = [A Zz ; Zz Zz];
%    B = [B_tt B_tz ; B_zt B_zz];

    A = A;

    B = B_tz*inv(B_zz)*B_zt - B_tt;

end % For main function, "cem2D_createABmat_lossless"

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
