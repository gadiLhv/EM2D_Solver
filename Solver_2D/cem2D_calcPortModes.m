function [kz,fc,Et,Ez,eIdxs] = cem2D_calcPortModes(meshData,meshProps,materialList,materialAssign,simProps,f_sim)
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

    [eIdxs,nEdges] = mesh2D_createEdgeIndexing(meshData);

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
        Ae = imr*Ee(eL,Det) - (er*k0^2)*Fe(eL,Det,b,c);

        % B matrix top left sub-matrix is still pretty simple
        Be_tt = imr*Fe(eL,Det,b,c);

        % Here it gets nasty
        eL = permute(eL,[2 3 1]);
        Det = permute(Det,[3 2 1]);
        b = permute(b,[2 3 1]);
        c = permute(c,[2 3 1]);

        % "f" elements defined in page 242, "Finite element method in electromagnetics", 1st edition.
        f = @(i,j) (b(i,:,:).*b(j,:,:) + c(i,:,:).*c(j,:,:));

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

        Be_zt = permute(Be_tz,[2 1 3]);

%        Be_zt = Be_tz;
%        Be_zt_11 = (imr*eL(1,:,:)./(Det*12)).*(f(1,2) - f(1,1));
%        Be_zt_12 = (imr*eL(2,:,:)./(Det*12)).*(f(1,3) - f(1,2));
%        Be_zt_13 = (imr*eL(3,:,:)./(Det*12)).*(f(1,1) - f(1,3));
%
%        Be_zt_21 = (imr*eL(1,:,:)./(Det*12)).*(f(2,2) - f(2,1));
%        Be_zt_22 = (imr*eL(2,:,:)./(Det*12)).*(f(2,3) - f(2,2));
%        Be_zt_23 = (imr*eL(3,:,:)./(Det*12)).*(f(2,1) - f(2,3));
%
%        Be_zt_31 = (imr*eL(1,:,:)./(Det*12)).*(f(3,2) - f(3,1));
%        Be_zt_32 = (imr*eL(2,:,:)./(Det*12)).*(f(3,3) - f(3,2));
%        Be_zt_33 = (imr*eL(3,:,:)./(Det*12)).*(f(3,1) - f(3,3));
%
%        Be_zt = [   Be_zt_11 Be_zt_12 Be_zt_13 ; ...
%                    Be_zt_21 Be_zt_22 Be_zt_23 ; ...
%                    Be_zt_31 Be_zt_32 Be_zt_33];

        % And finally, cross-pol elements:
        Be_zz = zeros([3 3 numel(Det)]);
        for i = 1:3
            for j = 1:3
                Be_zz(i,j,:) = imr*f(i,j)./(4*Det) - (k0^2)*er*(1 + (i == j)).*Det/12;
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

    % Find all edges existant in only one triangle, namely, boundary nodes
    uIdxs = 1:nEdges;
    howManyTris = sum(sum(bsxfun(@eq,permute(eIdxs,[1 3 2]),uIdxs),3),1);

    boundaryEdges = find(howManyTris(:) == 1);
    nonBoundaryEdges = uIdxs;
    nonBoundaryEdges(boundaryEdges) = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Debug: Find boundary nodes, to impose Dirichlet boundary conditions %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    [inWhichTri,whichEdge] = find(sum(eIdxs == permute(boundaryEdges,[3 2 1]),3));
%    figure;
%    for tIdx = unique(meshData.tnum).'
%        patch(  'faces',meshData.tria(meshData.tnum == tIdx,1:3),'vertices',meshData.vert, ...
%                'facecolor','none', ...
%                'edgecolor',[0,0,0]) ;
%    end
%
%    hold on;
%
%    tris = meshData.tria;
%    for beIdx = 1:numel(inWhichTri)
%        cTriNode = tris(inWhichTri(beIdx),whichEdge(beIdx));
%        nTriNode = tris(inWhichTri(beIdx),mod(whichEdge(beIdx),3) + 1);
%
%        cTriNode = meshData.vert(cTriNode,:);
%        nTriNode = meshData.vert(nTriNode,:);
%
%        plot([cTriNode(1) ; nTriNode(1)],[cTriNode(2) ; nTriNode(2)],'-r','linewidth',3);
%    end
%    hold off;
%    close(gcf);
    %%%%%%%%%%%%%
    % End Debug %
    %%%%%%%%%%%%%

%    A = removeBoundaryEdges(A,boundaryEdges);
%    B_tt = removeBoundaryEdges(B_tt,boundaryEdges);
%    B_tz = removeBoundaryEdges(B_tz,boundaryEdges);
%    B_zt = removeBoundaryEdges(B_zt,boundaryEdges);
%    B_zz = removeBoundaryEdges(B_zz,boundaryEdges);


    %%%%%%%%%%%%
    % Option 1 %
    %%%%%%%%%%%%

%    % On edges that are PEC (single sided edges), remove fields
%    % Option 1: Only Et
%    A = A;
%    B = B_tz*inv(B_zz)*B_zt - B_tt;
%
%%    A = removeBoundaryEdges(A,boundaryEdges);
%%    B = removeBoundaryEdges(B,boundaryEdges);
%
%
%    [Et,lambda] = eig(B\A,'vector');
%
%    Et_pad = zeros([nEdges size(Et,2)]);
%    Et_pad(nonBoundaryEdges,:) = Et;
%    Et = Et_pad;
%    clear('Et_pad');
%
%    kz2 = lambda;
%    kc2 = (2*pi*f_sim/c0)^2 - kz2;
%    fc = sqrt(kc2)*c0/(2*pi);
%
%    % Dump all imaginary fc (namely, negative wavenumbers)
%    fc_Th = 0.01;
%    binBadFc = abs(imag(fc)) > fc_Th;
%%    binBadFc = imag(fc) ~= 0;
%    kc2(binBadFc) = [];
%    kz2(binBadFc) = [];
%    Et(:,binBadFc) = [];
%    fc(binBadFc) = [];
%
%    fc = real(fc);
%
%    [fc,sortIdxs] = sort(real(fc));
%    kz2 = kz2(sortIdxs);
%    kc2 = kc2(sortIdxs);
%    Et = Et(:,sortIdxs);
%
%    kz = sqrt(kz2);
%    kz(imag(kz) > 0) = real(kz(imag(kz) > 0)) - 1i*imag(kz(imag(kz) > 0));
%
%    Et = Et./(kz(:).');
%
%    Ez = NaN;

    %%%%%%%%%%%%
    % Option 2 %
    %%%%%%%%%%%%

    Zz = zeros(size(A));
    A = [A Zz ; Zz Zz];
    B = [B_tt B_tz ; B_zt B_zz];

    [EtEz,lambda] = eig(B\A,'vector');

    Et = EtEz(1:(size(EtEz,1)/2),:);
    Ez = EtEz(((size(EtEz,1)/2)+1):end,:);

%    Et_pad = zeros([nEdges size(Et,2)]);
%    Et_pad(nonBoundaryEdges,:) = Et;
%    Et = Et_pad;
%    Ez_pad = zeros([nEdges size(Ez,2)]);
%    Ez_pad(nonBoundaryEdges,:) = Ez;
%    Ez = Ez_pad;
%    clear('Et_pad','Ez_pad');

    % Delete first half
    lambda(1:(size(Ez,1))) = [];
    Et(:,1:(size(Ez,1))) = [];
    Ez(:,1:(size(Ez,1))) = [];

    kz2 = -real(lambda);
    kc2 = (2*pi*f_sim/c0)^2 - kz2;

    badFc = kc2 < 0;
    kc2 = kc2(~badFc);
    Et = Et(:,~badFc);
    Ez = Ez(:,~badFc);

    fc = sqrt(kc2)*c0/(2*pi);

    % Dump all imaginary fc (namely, negative wavenumbers)
    fc_Th = 0.01;
    binBadFc = abs(imag(fc)) > fc_Th;
%    binBadFc = imag(fc) ~= 0;
    kc2(binBadFc) = [];
    kz2(binBadFc) = [];
    Et(:,binBadFc) = [];
    Ez(:,binBadFc) = [];
    fc(binBadFc) = [];

    fc = real(fc);

    [fc,sortIdxs] = sort(real(fc));
    kz2 = kz2(sortIdxs);
    kc2 = kc2(sortIdxs);
    Et = Et(:,sortIdxs);

    kz = sqrt(kz2);
    kz(imag(kz) > 0) = real(kz(imag(kz) > 0)) - 1i*imag(kz(imag(kz) > 0));

    Et = Et./(kz(:).');

    Ez = Ez/(-1i);

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
