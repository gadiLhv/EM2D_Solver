function [a,b,c,d,V] = cem3D_createInterpolantCoeffs1st(vert,tet)

    X = [vert(tet(:,1),1) vert(tet(:,2),1) vert(tet(:,3),1) vert(tet(:,4),1)]
    Y = [vert(tet(:,1),2) vert(tet(:,2),2) vert(tet(:,3),2) vert(tet(:,4),2)]
    Z = [vert(tet(:,1),3) vert(tet(:,2),3) vert(tet(:,3),3) vert(tet(:,4),3)]

    x = @(i) permute(X(:,i),[3 2 1]);
    y = @(i) permute(Y(:,i),[3 2 1]);
    z = @(i) permute(Z(:,i),[3 2 1]);
    
    A1 = [[x(2) x(3) x(4)] ; [y(2) y(3) y(4)] ; [z(2) z(3) z(4)]];
    A2 = [[x(1) x(3) x(4)] ; [y(1) y(3) y(4)] ; [z(1) z(3) z(4)]];
    A3 = [[x(1) x(2) x(4)] ; [y(1) y(2) y(4)] ; [z(1) z(2) z(4)]];
    A4 = [[x(1) x(2) x(3)] ; [y(1) y(2) y(3)] ; [z(1) z(2) z(3)]];
    
    I = ones([1 1 size(X,1)]);
    
    B1 = [[I  I  I] ; [y(2) y(3) y(4)] ; [z(2) z(3) z(4)]];
    B2 = [[I  I  I] ; [y(1) y(3) y(4)] ; [z(1) z(3) z(4)]];
    B3 = [[I  I  I] ; [y(1) y(2) y(4)] ; [z(1) z(2) z(4)]];
    B4 = [[I  I  I] ; [y(1) y(2) y(3)] ; [z(1) z(2) z(3)]];
    
    C1 = [[x(2) x(3) x(4)] ; [I  I  I] ; [z(2) z(3) z(4)]];
    C2 = [[x(1) x(3) x(4)] ; [I  I  I] ; [z(1) z(3) z(4)]];
    C3 = [[x(1) x(2) x(4)] ; [I  I  I] ; [z(1) z(2) z(4)]];
    C4 = [[x(1) x(2) x(3)] ; [I  I  I] ; [z(1) z(2) z(3)]];
    
    D1 = [[x(2) x(3) x(4)] ; [y(2) y(3) y(4)] ; [I  I  I]];
    D2 = [[x(1) x(3) x(4)] ; [y(1) y(3) y(4)] ; [I  I  I]];
    D3 = [[x(1) x(2) x(4)] ; [y(1) y(2) y(4)] ; [I  I  I]];
    D4 = [[x(1) x(2) x(3)] ; [y(1) y(2) y(3)] ; [I  I  I]];
    
    v = [[I I I I] ; [x(1) x(2) x(3) x(4)] ; [y(1) y(2) y(3) y(4)] ; [z(1) z(2) z(3) z(4)]];
    
    
    assembleMatrix = @(m1,m2,m3,m4) permute([computeDet3x3(m1) computeDet3x3(m2) computeDet3x3(m3) computeDet3x3(m4)],[3 2 1]);
    
    a = assembleMatrix(A1,A2,A3,A4);
    b = -assembleMatrix(B1,B2,B3,B4);
    c = -assembleMatrix(C1,C2,C3,C4);
    d = -assembleMatrix(D1,D2,D3,D4);
    
    V = computeMinorDet4x4(v,1,1) - computeMinorDet4x4(v,1,2) + computeMinorDet4x4(v,1,3) - computeMinorDet4x4(v,1,4);
    V = permute(V,[3 2 1]);
end


function Dij = computeMinorDet4x4(A,i,j)
    % Determine minor indices
    I = 1:4;
    J = 1:4;
    I(i) = [];
    J(j) = []; 
    Mij = A(I,J,:);
    
    Dij = computeDet3x3(Mij);
end


function Dij = computeDet3x3(Mij)
    Dij = Mij(1,1,:).*(Mij(2,2,:).*Mij(3,3,:) - Mij(2,3,:).*Mij(3,2,:)) + ...
          Mij(1,2,:).*(Mij(2,3,:).*Mij(3,1,:) - Mij(2,1,:).*Mij(3,3,:)) + ...
          Mij(1,3,:).*(Mij(2,1,:).*Mij(3,2,:) - Mij(2,2,:).*Mij(3,1,:));
end
