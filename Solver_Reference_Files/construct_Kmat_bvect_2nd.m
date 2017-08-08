function [K,b] = construct_Kmat_bvect_2nd(nodeList,elementList,segmentList,ax,ay,beta,gammas1,gammas2,qs)

%% Phrase the Ke (element) matrix, excluding the edges

% This will be a generalized calculation of all the coefficients. The index
% numbers for the vertices will be global index numbers. A matrix will be
% constructed for each object and added to the previous one

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

ke = zeros([3 3 size(ci,1)]);

for i = 1:3
    for j = 1:3
        ke(i,j,:) = ((4*dete).^(-1)).*(ax.*bi(:,i).*bi(:,j) + ay.*ci(:,i).*ci(:,j)) + (dete/12).*beta.*(1 + (i == j));
    end
end

%% Calculate boundaries (segments) as Ks matrices and bs vectors

segStartPos = nodeList(segmentList(:,2),2:3);
segEndPos = nodeList(segmentList(:,3),2:3);

% Calculate difference vectors
drs = segEndPos - segStartPos;
% Sizes of segments
ls = sqrt(sum(drs.^2,2));

% ks matrix elements
ks = zeros([2 2 size(ls,1)]);

for i = 1:2
    for j = 1:2
        ks(i,j,:) = gammas1.*ls/6 - gammas2./ls + (i==j)*(gammas1.*ls/6 + 2*gammas2./ls);
    end
end

% bs vector elements
bs = zeros([2 size(ls,1)]);



for i = 1:2
    bs(i,:) = qs.*ls/2;
end

%% Generalize matrix to K matrix and b vector

K = zeros([size(nodeList,1) size(nodeList,1)]);
b = zeros([size(nodeList,1) 1]);

% Add the elements from the local Ke matrices
for e = 1:size(ke,3)
    vert = [elementList(e,2) elementList(e,3) elementList(e,4)];
    K(vert,vert) = K(vert,vert) + ke(:,:,e);
end

% Add the elements from the local Ks matrics
for s = 1:size(ks,3)
    vert = [segmentList(s,2) segmentList(s,3)];
    K(vert,vert) = K(vert,vert) + ks(:,:,s);
end

% Add the elements from the local bs vectors
for s = 1:size(ks,3)
    vert = [segmentList(s,2) segmentList(s,3)];
    b(vert) = b(vert) + bs(:,s);
end