function [xMesh,yMesh,phiI] = thickenMeshByElement(nodeList,elementList,phi,range,nMesh)

xi = [nodeList(elementList(:,2),2) nodeList(elementList(:,3),2) nodeList(elementList(:,4),2)];
yi = [nodeList(elementList(:,2),3) nodeList(elementList(:,3),3) nodeList(elementList(:,4),3)];

ai = [  xi(:,2).*yi(:,3) - yi(:,2).*xi(:,3) ...
        xi(:,3).*yi(:,1) - yi(:,3).*xi(:,1) ...
        xi(:,1).*yi(:,2) - yi(:,1).*xi(:,2)];
bi = [  yi(:,2) - yi(:,3) ...
        yi(:,3) - yi(:,1) ...
        yi(:,1) - yi(:,2)];
ci = [  xi(:,3) - xi(:,2) ...
        xi(:,1) - xi(:,3) ...
        xi(:,2) - xi(:,1)];

dete = 0.5*(bi(:,1).*ci(:,2)- bi(:,2).*ci(:,1));

[xMesh,yMesh] = meshgrid(linspace(range(1),range(2),nMesh(1)),...
                        linspace(range(3),range(4),nMesh(2)));

% phiI = zeros(size(xMesh));
phiI = ones(size(xMesh))*NaN;

for e = 1:size(elementList,1)
    vert = nodeList(elementList(e,2:4)',2:3);
    vertVals = phi(elementList(e,2:4)');
    [inPoly,onPoly] = inpolygon(xMesh,yMesh,vert(:,1),vert(:,2));
    binArr = inPoly | onPoly;
    
    xList = xMesh(binArr);
    yList = yMesh(binArr);
    
    phii = zeros([numel(xList) 1]);
    
    for j = 1:3
        phii = phii + ((2*dete(e)).^(-1)).*((ai(e,j) + ...
                    bi(e,j)*xList + ...
                    ci(e,j)*yList)*vertVals(j));
    end
    
    phiI(binArr) = phii;
    
end