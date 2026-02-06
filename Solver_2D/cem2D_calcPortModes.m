function [lambda,Et,Ez,eIdxs] = cem2D_calcPortModes(meshData,meshProps,materialList,materialAssign,simProps,f_sim)
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
    nVerts = size(vert_m,1);

    [eIdxs,nEdges,edge2vert] = mesh2D_createEdgeIndexing(meshData);

    A_tt = sparse(nEdges, nEdges);
    B_tt = sparse(nEdges, nEdges);
    B_tz = sparse(nEdges, nVerts);
    B_zt = sparse(nVerts, nEdges);
    B_zz = sparse(nVerts, nVerts);

    % Find all edges existant in only one triangle, namely, boundary nodes
    uIdxs = 1:nEdges;
    howManyTris = sum(sum(bsxfun(@eq,permute(eIdxs,[1 3 2]),uIdxs),3),1);

    boundaryEdges = find(howManyTris(:) == 1);
    nonBoundaryEdges = uIdxs;
    nonBoundaryEdges(boundaryEdges) = [];

    % Find all boundary vertices
    boundaryVerts = edge2vert(boundaryEdges,:);
    boundaryVerts = unique(boundaryVerts(:));

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
        eL = cem2D_calcEdgeLengths(vert_m,triTriplets);

%        % Create interpolants for elements
        [a,b,c,Det] = cem2D_createInterpolantCoeffs1st(vert_m,triTriplets);

        gradN_x = b./(2*Det);
        gradN_y = c./(2*Det);

        curlN = calc_curlN(eL,Det);

        % % A matrix is simple enough
        % Ae = imr*Ee(eL,Det) - (er*k0^2)*Fe(eL,Det,b,c);

        % % B matrix top left sub-matrix is still pretty simple
        % Be_tt = imr*Fe(eL,Det,b,c);

        % Modification on the LLM variation. Calculates same matrix (hopefully)
        % but more explicitly written
        A_tt_e = zeros(3,3,numel(Det));
        B_tt_e = zeros(3,3,numel(Det));

        for i = 1:3
            for j = 1:3
                A_tt_e(i,j,:) = imr*curlN(:,i).*curlN(:,j).*Det ...
                                 - (k0^2)*er*edgeMass(i,j,Det);

                B_tt_e(i,j,:) = imr*edgeMass(i,j,Det);
            end
        end

        B_tz_e = zeros(3,3,numel(Det));
        for i = 1:3
            for j = 1:3
                B_tz_e(i,j,:) = imr*...
                  (eL(:,i)./(2*Det)).*...
                  ( (b(:,mod(i,3)+1)-b(:,i)).*gradN_x(:,j) + ...
                    (c(:,mod(i,3)+1)-c(:,i)).*gradN_y(:,j)).*Det/6;
            end
        end
        B_zt_e = permute(B_tz_e,[2 1 3]);

        % % And finally, cross-pol elements:
        % Be_zz = zeros([3 3 numel(Det)]);
        % for i = 1:3
        %     for j = 1:3
        %         Be_zz(i,j,:) = imr*f(i,j)./(4*Det) - (k0^2)*er*(1 + (i == j)).*Det/12;
        %     end
        % end
        B_zz_e = zeros(3,3,numel(Det));
        for i = 1:3
          for j = 1:3
            B_zz_e(i,j,:) = ...
              imr*(gradN_x(:,i).*gradN_x(:,j) + gradN_y(:,i).*gradN_y(:,j)).*Det ...
              -(k0^2)*er*nodalMass(i,j,Det);
          end
        end


        % Iteratively add sub-matrices
        for e = 1:size(A_tt_e,3)
            ceIdxs = edgeTripliets(e,:);
            cvIdxs = meshData.tria(e,:);

            A_tt(ceIdxs,ceIdxs) = A_tt(ceIdxs,ceIdxs) + A_tt_e(:,:,e);
            B_tt(ceIdxs,ceIdxs) = B_tt(ceIdxs,ceIdxs) + B_tt_e(:,:,e);
            B_tz(ceIdxs,cvIdxs) = B_tz(ceIdxs,cvIdxs) + B_tz_e(:,:,e);
            B_zt(cvIdxs,ceIdxs) = B_zt(cvIdxs,ceIdxs) + B_zt_e(:,:,e);
            B_zz(cvIdxs,cvIdxs) = B_zz(cvIdxs,cvIdxs) + B_zz_e(:,:,e);
        end

    end % Per-face for loop

    %%%%%%%%%%%%%
    % End Debug %
    %%%%%%%%%%%%%

%    A = removeBoundaryEdges(A,boundaryEdges);
%    B_tt = removeBoundaryEdges(B_tt,boundaryEdges);
%    B_tz = removeBoundaryEdges(B_tz,boundaryEdges);
%    B_zt = removeBoundaryEdges(B_zt,boundaryEdges);
%    B_zz = removeBoundaryEdges(B_zz,boundaryEdges);

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
%
%
%    sHdl = scatter(meshData.vert(boundaryVerts,1), meshData.vert(boundaryVerts,2),20,'filled');
%
%    hold off;
%    close(gcf);
    %%%%%%%%%%%%%
    % End Debug %
    %%%%%%%%%%%%%

    A_tz = zeros(size(B_tz));
    A_zt = zeros(size(B_zt));
    A_zz = zeros(size(B_zz));
    B = [B_tt B_tz ; B_zt B_zz];
    A = [A_tt A_tz ; A_zt A_zz];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Option 1: Buld matrices with boundary conditions inherent %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    binVertsEdges = true([(nVerts + nEdges) 1]);
    binVertsEdges(boundaryEdges) = false;
    binVertsEdges(nEdges + boundaryVerts) = false;

    A = removeBoundaryEdges(A,~binVertsEdges);
    B = removeBoundaryEdges(B,~binVertsEdges);

    % End of options. Solve Eigenvalue problem
    [EtEz_part,lambda] = eig(A,B,'vector');

    EtEz = zeros([(nVerts + nEdges) size(EtEz_part,2)]);


    EtEz(binVertsEdges,:) = EtEz_part;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Option 2: Build matrices, unite them, then add BC %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%    A = inv(B)*A;
%    B = eye(size(A,1));
%
%    A(boundaryEdges + (boundaryEdges - 1)*size(A_tt,1)) = 1;
%    A(boundaryEdges + (boundaryVerts - 1)*size(B_tz,1)) = 1;
%    A(boundaryVerts + (boundaryEdges - 1)*size(B_zt,1)) = 1;
%    A(boundaryVerts + (boundaryVerts - 1)*size(B_zz,1)) = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Option 3: Impose Dirichlet B.C. explicitly with larger matrix %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    A_tz = zeros(size(B_tz));
%    A_zt = zeros(size(B_zt));
%    A_zz = zeros(size(B_zz));
%
%    A = [A_tt A_tz ; A_zt A_zz];
%    B = [B_tt B_tz ; B_zt B_zz];
%
%    Zz = zeros([1 1]*(nEdges + nVerts));
%    Ad = zeros([(nEdges + nVerts) 1]);
%    Ad(boundaryEdges) = 1;
%    Ad(boundaryVerts + nEdges) = 1;
%    Ad = diag(Ad);
%
%    % Enforce on boundary edges and vertices
%    Afull = [A Zz ; Zz Ad];
%    Bfull = [B Zz ; Zz Zz];
%
%    % End of options. Solve Eigenvalue problem
%    [EtEz,lambda] = eig(Afull,Bfull,'vector');
%
%    EtEz = EtEz(1:(nEdges + nVerts),1:(nEdges + nVerts));
%    lambda = lambda(1:(nEdges + nVerts));

    %%%%%%%%%%%%%%%
    % End options %
    %%%%%%%%%%%%%%%

%    EtEz = EtEz.';
    Et = EtEz(1:nEdges,:);
    Ez = EtEz((nEdges + 1):end,:);

    k02 = (2*pi*f_sim/c0)^2;
    kz2 = lambda;
    kc2 = k02 - kz2;

%    badFc_th = 0.001;
%    badFc = real(kc2) < 0;
%    badFc = (real(kc2) < 0) | (abs(imag(kc2)) ~= 0);
%    badFc = (real(kc2) < 0) | (abs(imag(kc2)) >= badFc_th) | (abs(real(kc2)) <= badFc_th);

%    lambda = lambda(~badFc);
%    kz2 = kz2(~badFc);
%    kc2 = kc2(~badFc);
%    Et = Et(:,~badFc);
%    Ez = Ez(:,~badFc);

%    binTE = binTE(~badFc);

%    fc = sqrt(kc2)*c0/(2*pi);

%    [fc,sortIdxs] = sort(real(fc));
%    kz2 = kz2(sortIdxs);
%    kc2 = kc2(sortIdxs);
%    Et = Et(:,sortIdxs);
%    binTE = binTE(sortIdxs);

    kz = sqrt(kz2);
    kz(imag(kz) > 0) = real(kz(imag(kz) > 0)) - 1i*imag(kz(imag(kz) > 0));

    Et = Et./(kz(:).');
    Ez = Ez/(-1i);

%    fc_TE = fc(binTE);
%    fc_TM = fc(~binTE);

end % For main function, "cem2D_createABmat_lossless"

function cE = Ee(edgeL,detE)

    edgeL = permute(edgeL,[2 3 1]);
    detE = permute(detE,[3 2 1]);

    % Calculate the Eij matrix per element
    cE = bsxfun(@times,bsxfun(@times,edgeL,permute(edgeL,[2 1 3])),1./detE);

end

function curlN = calc_curlN(Le,Det)
    curlN = Le./(2*Det);
end

function M = edgeMass(i,j,Det)
  if i == j
    M = Det/6;
  else
    M = Det/12;
  end
end

function M = nodalMass(i,j,Det)
  if i == j
    M = Det/6;
  else
    M = Det/12;
  end
end

function [gradN,b,c,Det] = calc_gradN(verts)
    % Geometry
    x = verts(:,1);
    y = verts(:,2);

    Det = abs(det([1 x(1) y(1);
                   1 x(2) y(2);
                   1 x(3) y(3)])) / 2;

    % Gradients of nodal basis
    b = [y(2)-y(3);
         y(3)-y(1);
         y(1)-y(2)];

    c = [x(3)-x(2);
         x(1)-x(3);
         x(2)-x(1)];

    gradN = [b c] / (2*Det);   % 3กั2
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

function A = removeBoundaryEdges(A,i)
    A(:,i) = [];
    A(i,:) = [];
end

