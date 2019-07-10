function [Hz,k_xi] = cem1D_calcTEmodes(segVerts,meshData,meshProps,segMaterials,materialList,simProps)
% [Hz,k_xi] = cem1D_calcTEmodes(segVerts,meshData,meshProps,segMaterials,materialList,simProps)
% Calculates the 1D modes, with Z-infinite approximation, for the TE case (Hz).
% Inputs:
% 1. segVerts - Segment vertices index pairs, as they appear in meshData.vert
% 2. meshData - Mesh data structure passed on from the mesh generation functions
% 3. meshProps - Properties passed on to the mesher.
% 4. segMaterials - Cell list indicating the material name of each segment.
% 5. materialList - Cell list that is filled with material list strucures.
% 6. simProps - Simulation properties strucures
%
% Outputs:
% 1. Hz - [N_segments,N_eigenmodes] complex array. Each column denotes an eigenmode.
% 2. k_xi - N_eigenmode element vector. Denoting the tangential wave number.
    
    % Extract metallic segments out of the equations.
    % The fields in these areas are zero, by identity. 
    getMaterials = @(cCell) cem2D_getMaterialPropsFromName(cCell,materialList);
    matList = cellfun(getMaterials,segMaterials);
    isNormalMaterial = @(cStruct) strcmp(cStruct.type,'normal');
    binIsNormal = arrayfun(isNormalMaterial,matList).';
    
    normalSegIdxs = find(binIsNormal);
    
    % Currently exceprt only normal materials
    origSegVerts = segVerts;
    segVerts = segVerts(normalSegIdxs,:);
    segMaterials = segMaterials(normalSegIdxs);
    
    % Generate mapping to a new, shorter, vertex list
    [uniqueSegVerts,~,segVertsMapping] = unique(segVerts(:));
    segVertsMapping = reshape(segVertsMapping,[],2);
    
    pt1 = meshData.vert(segVerts(:,1),:);
    pt2 = meshData.vert(segVerts(:,2),:);
    
    % Initialize eigenvalue matrices
    A = zeros(size(segVerts,1) + 1);
    B = A;
    
    l = sqrt(sum(abs(pt2 - pt1).^2,2));
    
    % Convert l to MKS
    l = units(simProps.lengthUnits,'m',l);
    % Choose only normal materials for the mode calculations.
    % TODO 
    % 1. Convert "weakest material, to give dielectric constant to relevant metals??
    % 2. Convert TE\TM mode builders to exclude fully metallic segments.
    
    
    
    for segIdx = 1:size(segVerts,1)
        cMaterial = cem2D_getMaterialPropsFromName(segMaterials{segIdx},materialList);
        matType = cMaterial.type;
        er = cMaterial.er;
        mr = cMaterial.mr;
        
        a = -(1/er);
        b = mr;
        
        le = l(segIdx);
        
        if strcmp(matType,'normal');
            % Mark metals, as they require 
            % alpha and beta coefficients
            Ae_BC = [];
            Be_BC = [];
            % Currently only handles metals as pure metal
        else
            Ae_BC = ones([2 1])*a/le;
            Be_BC = -ones([2 1])*b*le/3;
        end
        
        
        % Initialize sub-element matrices
        Ae = (1/le)*[[a -a] ; [-a a]];
        Be = -b*le*[[1/3 1/6] ; [1/6 1/3]];
        
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
   
    [Unormal,eigVal] = eig(B\A);
    k_xi = sqrt(diag(eigVal));
        
    % Sort by eigenvalue order
    [k_xi,sortIdxs] = sort(k_xi);
    Unormal = Unormal(:,sortIdxs);
        
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
    
    Hz = U;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Debug: Draw first 3 modes %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    figure('position',[257     75   1138    749]);
%    c0 = physical_constant('speed of light in vacuum');
%    e0 = physical_constant('electric constant');
%    m0 = physical_constant('mag. constant');
%    f0 = units(simProps.freqUnits,'Hz',simProps.fSim);
%    k0 = 2*pi*f0/c0;
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
%        ylabel('|H_z| [A/m]','fontsize',14);
%        kz = sqrt(k0^2 - k_xi(modeIdx).^2);
%        f_cutoff = c0*k_xi(modeIdx)/(2*pi);
%        er_eff = c0*kz./(2*pi*f0);
%        title(sprintf('Mode #%d ; f_c = %.2f GHz',modeIdx,f_cutoff/1e9),'fontsize',14);
%        grid on;
%    end
%    
%    close(gcf);
%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End Debug: Draw first 3 modes %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
