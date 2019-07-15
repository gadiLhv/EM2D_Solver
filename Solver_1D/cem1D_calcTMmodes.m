function [Ez,k_xi,segVerts] = cem1D_calcTMmodes(segVerts,meshData,meshProps,segMaterials,materialList,simProps)
% [Ez,k_xi] = cem1D_calcTMmodes(segVerts,meshData,meshProps,segMaterials,materialList,simProps)
% Calculates the 1D modes, with Z-infinite approximation, for the TM case (Ez).
% Inputs:
% 1. segVerts - Segment vertices index pairs, as they appear in meshData.vert
% 2. meshData - Mesh data structure passed on from the mesh generation functions
% 3. meshProps - Properties passed on to the mesher.
% 4. segMaterials - Cell list indicating the material name of each segment.
% 5. materialList - Cell list that is filled with material list strucures.
% 6. simProps - Simulation properties strucures
%
% Outputs:
% 1. Ez - [N_segments,N_eigenmodes] complex array. Each column denotes an eigenmode.
% 2. k_xi - N_eigenmode element vector. Denoting the tangential wave number.

    % Extract metallic segments out of the equations.
    % The fields in these areas are zero, by identity. 
    getMaterials = @(cCell) cem2D_getMaterialPropsFromName(cCell,materialList);
    matList = cellfun(getMaterials,segMaterials);
    isNormalMaterial = @(cStruct) strcmp(cStruct.type,'normal');
    binIsNormal = arrayfun(isNormalMaterial,matList).';
    
    % Mark segment indexes
    normalSegIdxs = find(binIsNormal);
    
    % For the TM case, also mark the metallic nodes
    metalSegIdxs = find(~binIsNormal);
    metalSegVerts = segVerts(metalSegIdxs,:);
    
    % Currently exceprt only normal materials
    origSegVerts = segVerts;
    segVerts = segVerts(normalSegIdxs,:);
    segMaterials = segMaterials(normalSegIdxs);
    
    
    % Generate mapping to a new, shorter, vertex list
    [uniqueSegVerts,~,segVertsMapping] = unique(segVerts(:));
    segVertsMapping = reshape(segVertsMapping,[],2);
    
    % Same for metal segments
    [uniqueMetalVerts,~,metalVertsMapping] = unique(metalSegVerts(:));
    metalVertsMapping = reshape(metalVertsMapping,[],2);
    
    pt1 = meshData.vert(segVerts(:,1),:);
    pt2 = meshData.vert(segVerts(:,2),:);
    
    % Initialize eigenvalue matrices
    A = zeros(size(segVerts,1) + 1);
    B = A;
    
    l = sqrt(sum(abs(pt2 - pt1).^2,2));
    l = units(simProps.lengthUnits,'m',l);
    
    % Same for metals
    pt1_metals = meshData.vert(origSegVerts(metalSegIdxs,1),:);
    pt2_metals = meshData.vert(origSegVerts(metalSegIdxs,2),:);
    l_metals = sqrt(sum(abs(pt2_metals - pt1_metals).^2,2));
    l_metals = units(simProps.lengthUnits,'m',l_metals);
    
    % Choose only normal materials for the mode calculations.
    % TODO 
    % 1. Convert "weakest material, to give dielectric constant to relevant metals??
    % 2. Convert TE\TM mode builders to exclude fully metallic segments.
    
    e0 = physical_constant('electric constant');
    m0 = physical_constant('mag. constant');
    c0 = physical_constant('speed of light in vacuum');
    
    k0 = 2*pi*units(simProps.freqUnits,'Hz',simProps.fSim)/c0;
    
    for segIdx = 1:size(segVerts,1)
        cMaterial = cem2D_getMaterialPropsFromName(segMaterials{segIdx},materialList);
        matType = cMaterial.type;
        er = cMaterial.er;
        mr = cMaterial.mr;
        
        a = -(1/mr);
        b = er*k0^2;   
        g = 1/mr;
        
        le = l(segIdx);
        
        if strcmp(matType,'normal');
            % Mark metals, as they require 
            % alpha and beta coefficients
            Ae_BC = [];
            Be_BC = [];
            % Currently only handles metals as pure metal
        else
%            Ae_BC = ones([2 1])*a/le;
%            Be_BC = -ones([2 1])*b*le/3;
            Ae_BC = ones([2 1])*a/le - ones([2 1])*b*le/3;
            Be_BC = -ones([2 1])*g*le/3;
        end
        
        % Initialize sub-element matrices
        Ae = (1/le)*[[a -a] ; [-a a]] + b*le*[[1/3 1/6] ; [1/6 1/3]];
        Be = g*le*[[1/3 1/6] ; [1/6 1/3]];
        
        % Here "hard" B.C. of zero electric TANGENTIAL electric field
        % must be suited inside.
        
        segRange = segVertsMapping(segIdx,:);
        A(segRange,segRange) = A(segRange,segRange) + Ae;
        B(segRange,segRange) = B(segRange,segRange) + Be;
        
        % Add neumann boundary conditions of a PEC (dHz/dx = 0)
        if ~isempty(Ae_BC)
            for idx = 1:2
                A(segRange(idx),segRange(idx)) = A(segRange(idx),segRange(idx)) + Ae_BC(idx);
                B(segRange(idx),segRange(idx)) = B(segRange(idx),segRange(idx)) + Be_BC(idx);
            end
        end
        
        
    end

    % Force Ez = 0 at metals by removing the relevant rows and columns
    % 1. Find the metal vertices in the "unique metal vertices"
    [idxInNormals,~] = find(bsxfun(@eq,uniqueSegVerts(:),metalSegVerts(:).'));
    returnIdxs = (1:numel(uniqueSegVerts)).';
    returnIdxs(idxInNormals) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Debug: Show points to be removed %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    figure;
%    plot(pt1(:,1),pt1(:,2),'.b','markersize',25);
%    hold on;
%    plot(pt2(:,1),pt2(:,2),'.b','markersize',25);
%    
%    segCoors = meshData.vert(uniqueSegVerts,:);
%    
%    plot(pt1_metals(:,1),pt1_metals(:,2),'.g','markersize',17);
%    plot(pt2_metals(:,1),pt2_metals(:,2),'.g','markersize',17);
%    
%    plot(segCoors(idxInNormals,1),segCoors(idxInNormals,2),'or','markersize',17);
%    hold off;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End Debug: Show points to be removed %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A_orig = A;
    B_orig = B;
    
    A(idxInNormals,:) = [];
    B(idxInNormals,:) = [];
    A(:,idxInNormals) = [];
    B(:,idxInNormals) = [];
   
    [Unormal,eigVal] = eig(B\A);
    k_zeta = sqrt(diag(eigVal));
    k_xi = sqrt(k0^2 - k_zeta.^2);
    
    if(max(abs(imag(k_xi(:)))) > 1e-6)
        warning('Cutoff frequency is imaginary up to 1e-6.');
    end
    
    [k_xi,sortIdxs] = sort(real(k_xi));
    Unormal = Unormal(:,sortIdxs);
    
    % Pad the Unormal back
    Unormal_padded = zeros([size(A_orig,1) size(A,2)]);
    Unormal_padded(returnIdxs,:) = Unormal;
    Unormal = Unormal_padded;
    
    % Now this needs to be mapped into the original segment
    % structure. 
    % Dim. #1: Segment 
    % Dim. #2: Vertices of segments
    # Dim. #3: Eigenvalue
    # This is only the normal material list. 
    Unormal = permute(Unormal,[1 3 2]);
    Unormal = reshape(Unormal(segVertsMapping,:,:),size(segVerts,1),2,size(Unormal,3));
    
    % TODO: Finish mapping to original indexing.
    U = zeros([size(origSegVerts,1) 2 numel(k_xi)]);
    U(normalSegIdxs,:,:) = Unormal;
    segVerts = origSegVerts;
    
    Ez = U;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Debug: Draw first 3 modes %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    figure('position',[257     75   1138    749]);
%    
%%    c0 = physical_constant('speed of light in vacuum');
%%    f0 = units(simProps.freqUnits,'Hz',simProps.fSim);
%    
%    nSubX = 3;
%    nSubY = 2;
%    
%    % Find first and last point
%    pt1 = meshData.vert(segVerts(:,1),:);
%    pt2 = meshData.vert(segVerts(:,2),:);
%    distMap = sqrt(bsxfun(@minus,pt1(:,1),pt2(:,1).').^2 + bsxfun(@minus,pt1(:,2),pt2(:,2).').^2);
%    [L,maxDistIdx] = max(distMap(:));
%    [pt1idx,pt2idx] = ind2sub(size(segVerts,1)*[1 1],maxDistIdx);
%    % Find normalized direction (for projection)
%    lineDir = (pt2(pt2idx,:) - pt1(pt1idx,:));
%    lineDir = lineDir./sqrt(sum(lineDir.^2,2));
%    % Project the points on this line
%    xi1 = sum(bsxfun(@times,pt1 - pt1(pt1idx,:),lineDir),2);
%    xi2 = sum(bsxfun(@times,pt2 - pt1(pt1idx,:),lineDir),2);
%    
%    for modeIdx = 1:size(U,3)
%        subplot(nSubX,nSubY,modeIdx);
%        for segIdx = 1:size(xi1,1)
%            cU = U(segIdx,:,modeIdx).';    
%            plot([xi1(segIdx) ; xi2(segIdx)],abs(cU),'-b','linewidth',3);
%            hold on;
%        end
%        hold off;
%        set(gca,'fontsize',14);
%        xlabel('\xi [mm]','fontsize',14);
%        ylabel('|E_z| [V/m]','fontsize',14);
%        kz = sqrt(k0^2 - k_xi(modeIdx).^2);
%        f_cutoff = c0*k_xi(modeIdx)/(2*pi);
%%        er_eff = c0*kz./(2*pi*f0);
%        title(sprintf('Mode #%d ; f_c = %.2f GHz',modeIdx,f_cutoff/1e9),'fontsize',14);
%        grid on;
%    end
%    
%    close(gcf);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End Debug: Draw first 3 modes %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
