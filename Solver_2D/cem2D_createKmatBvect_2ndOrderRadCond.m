function [K,b] = cem2D_createKmatBvect_2ndOrderRadCond(meshData,meshProps,radEdgeIdxs,materialList,materialAssign,simProps,f_sim)

    % Inputs: 
    % 1. meshData - Structure with all of the mesh data
    % 2. meshProps - General mesh properties
    % 3. radEdgeIdxs - Edge numbers in the original faces - meshData.edge(:,edgeIdx) 
    % 4. materialList - Cell array with all the material data structures.
    % 5. materialAssign - Material assignements per face (object).
    % 6. simProps - General simulation properties
    % 7. f_sim - Frequency (in simulation units) of current solution

    distTH = meshProps.stitchingTolerance;

    c0 = physical_constant('speed of light in vacuum');
    e0 = physical_constant('electric constant');
    m0 = physical_constant('mag. constant');

    % Convert light speed to correct units
    f_sim = units(simProps.freqUnits,'Hz',f_sim);
    % Calculate wave number
    k0 = 2*pi*f_sim/c0;

    % The edges are referred to from meshData.etri
    edgePairs = meshData.etri(radEdgeIdxs,:);

    % Extract points
    p1 = units(simProps.lengthUnits,'m',meshData.vert(edgePairs(:,1),:));
    p2 = units(simProps.lengthUnits,'m',meshData.vert(edgePairs(:,2),:));

    % Extract segment materials, one by one
    er_segs = zeros([size(p1,1) 1]);
    mr_segs = er_segs;
    for segIdx = 1:size(p1,1)
        cSegMaterial = misc_segmentWeakestMaterial(p1(segIdx,:),p2(segIdx,:),meshData, materialList, materialAssign, meshProps);
        cMaterialProps = cem2D_getMaterialPropsFromName(cSegMaterial,materialList);    
        er_segs(segIdx) = cMaterialProps.er;
        mr_segs(segIdx) = cMaterialProps.mr;
    end

    nVerts = size(meshData.vert,1);

    % Initial FEM matrix and source vector
    K = zeros([1 1]*nVerts);
    b = zeros([nVerts 1]);

    % Calculate absoulte distance from axis origin.
    rho1 = meshData.vert(edgePairs(:,1),:);
    rho2 = meshData.vert(edgePairs(:,2),:);
    
    % Calculate Kappa coefficients for 'gamma1\2'
    kappa = 1./[sqrt(sum(rho1.^2,2)) sqrt(sum(rho2.^2,2))];

%    % Edit attempt: Curvature is always zero
%    kappa = zeros([size(edgePairs,1) 1]);

    k = k0*sqrt(er_segs.*mr_segs);
    % Calculate average gamma values
    gamma1 = 1i*k + kappa/2 - 1i*(kappa.^2)./(8*(1i*kappa - k0));
    gamma2 = -1i./(2*(1i*kappa - k));
    
    gamma1 = 0.5*sum(gamma1,2);
    gamma2 = 0.5*sum(gamma2,2);    
            
    % Calculate segment lengths
    l_s = sqrt(sum((rho2 - rho1).^2,2));

    % Construct sub matrices
    Ks = cat(3,[(gamma1.*l_s/3 + gamma2./l_s) (gamma1.*l_s/6 - gamma2./l_s)],...
               [(gamma1.*l_s/6 - gamma2./l_s) (gamma1.*l_s/3 + gamma2./l_s)]);

    % Transform both matrices to matrix batches
    Ks = permute(Ks,[3 2 1]);

    % Add the elements from the local Ks matrics
    for edgeIdx = 1:size(Ks,3)
        vertIdx = edgePairs(edgeIdx,:);
        K(vertIdx,vertIdx) = K(vertIdx,vertIdx) + Ks(:,:,edgeIdx);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculation for later - Plane-Wave port %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    % Segment centers
%    segCenters = (p1 + p2)./2;
%    % The angle is allready the mean value
%    segAngles = atan2(segCenters(:,2),segCenters(:,1));
%    
%    % Edit attempt - Angle of segment
%    segAngles = atan2(p2(:,2) - p1(:,2),p2(:,1) - p1(:,1));
%    
%    
%    q = exp(1i*k0*Rabs*cos(segAngles)-phi0)).*(...
%        -1i*k0*cos(segAngles(outerBin)-phi0) + 1i*k0 + kappa/2 - 1i*kappa./(8*(1i*kappa -  k0)) + ...
%        (1i./(2*(1i*kappa - k0))).*(k0.^2).*cos(segAngles - phi0).^2);


end

function binAns = inTH(a,b,TH)
    binAns = abs(a - b) < TH;
end