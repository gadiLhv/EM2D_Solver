function [fc_TE,fc_TM,E_TE,H_TM,eIdxs] = cem2D_calcModesByCutoff_TETM(meshData,meshProps,materialList,materialAssign,simProps)
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
%    f_sim = units(simProps.freqUnits,'Hz',f_sim);
    % Calculate wave number
%    k0 = 2*pi*f_sim/c0;

    vert_m = units(simProps.lengthUnits,'m',meshData.vert);

    [eIdxs,nEdges] = mesh2D_createEdgeIndexing(meshData);

    % Initial FEM matrix and source vector
    A_E = zeros([1 1]*nEdges);
    A_H = A_E;
    B_E = A_E;
    B_H = A_H;

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

        % Tangential component matrices

        % A matrix is simple enough
        Ae_E = (1/mr)*Ee(eL,Det);
        Ae_H = (1/er)*Ee(eL,Det);

        % B matrix top left sub-matrix is still pretty simple
        Be_E = er*Fe(eL,Det,b,c);
        Be_H = mr*Fe(eL,Det,b,c);

        % Here it gets nasty
        eL = permute(eL,[2 3 1]);
        Det = permute(Det,[3 2 1]);
        b = permute(b,[2 3 1]);
        c = permute(c,[2 3 1]);

        % "f" elements defined in page 242, "Finite element method in electromagnetics", 1st edition.
        f = @(i,j) (b(i,:,:).*b(j,:,:) + c(i,:,:).*c(j,:,:));

        % Iteratively add sub-matrices
        for e = 1:size(Ae_E,3)
            ceIdxs = edgeTripliets(e,:);
            A_E(ceIdxs,ceIdxs) = A_E(ceIdxs,ceIdxs) + Ae_E(:,:,e);
            A_H(ceIdxs,ceIdxs) = A_H(ceIdxs,ceIdxs) + Ae_H(:,:,e);
            B_E(ceIdxs,ceIdxs) = B_E(ceIdxs,ceIdxs) + Be_E(:,:,e);
            B_H(ceIdxs,ceIdxs) = B_H(ceIdxs,ceIdxs) + Be_H(:,:,e);
        end

    end % Per-face for loop

    % Find all edges existant in only one triangle, namely, boundary nodes
    uIdxs = 1:nEdges;
    howManyTris = sum(sum(bsxfun(@eq,permute(eIdxs,[1 3 2]),uIdxs),3),1);

    boundaryEdges = find(howManyTris(:) == 1);
    nonBoundaryEdges = uIdxs;
    nonBoundaryEdges(boundaryEdges) = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove boundary nodes
    %%%%%%%%%%%%%%%%%%%%%%%%%

%    A_E = removeBoundaryEdges(A_E,boundaryEdges);
%    B_E = removeBoundaryEdges(B_E,boundaryEdges);
%    A_H = removeBoundaryEdges(A_H,boundaryEdges);
%    B_H = removeBoundaryEdges(B_H,boundaryEdges);

    %%%%%%%%%%%%
    % Option 2 %
    %%%%%%%%%%%%

    Zz = zeros(size(A_E));
    A = [A_E Zz ; Zz A_H];
    B = [B_E Zz ; Zz B_H];

    [EteHtm,lambda] = eig(B\A,'vector');

    Ete = EteHtm(1:(size(EteHtm,1)/2),:);
    Htm = EteHtm(((size(EteHtm,1)/2)+1):end,:);

%    Ete_pad = zeros([nEdges size(Ete,2)]);
%    Ete_pad(nonBoundaryEdges,:) = Ete;
%    Ete = Ete_pad;
%    Htm_pad = zeros([nEdges size(Htm,2)]);
%    Htm_pad(nonBoundaryEdges,:) = Htm;
%    Htm = Htm_pad;
%    clear('Ete_pad','Htm_pad');

    % Binary vector that denotes if solution is TE or TM
    binTE = (1:size(Ete,2)).' <= size(Ete,2)/2;

    kc2 = lambda;

    badFc_th = 0.1;
    badFc = (real(kc2) < 0) | (abs(imag(kc2)) ~= 0);
%    badFc = real(kc2) < 0;
    kc2 = kc2(~badFc);
    Ete = Ete(:,~badFc);
    Htm = Htm(:,~badFc);
    binTE = binTE(~badFc);

    fc = sqrt(real(kc2))*c0/(2*pi);

    [fc,sortIdxs] = sort(fc);
    kc2 = kc2(sortIdxs);
    Ete = Ete(:,sortIdxs);
    Htm = Htm(:,sortIdxs);
    binTE = binTE(sortIdxs);

    E_TE = Ete(:,binTE);
    H_TM = Htm(:,~binTE);
    fc_TE = fc(binTE);
    fc_TM = fc(~binTE);

end % For main function, "cem2D_createABmat_lossless"

function cE = Ee(edgeL,detE)

    edgeL = permute(edgeL,[2 3 1]);
    detE = permute(detE,[3 2 1]);

    % Calculate the Eij matrix per element
    cE = bsxfun(@times,bsxfun(@times,edgeL,permute(edgeL,[2 1 3])),1./detE);

end

function A = removeBoundaryEdges(A,i)
    A(:,i) = [];
    A(i,:) = [];
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
