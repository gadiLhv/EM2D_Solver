function [nodeList,elementList,segmentList,N] = meshRect_wHobstacle(rectDims,maxD,obsDims,obsCenter)

% Inputs:
% rectDims - 2 element vector denoting the rectangle 
%            dimensions in meters [h,w]
% N -        2 element vector denoting the number of points
%            to dismember the side to [Nh,Nw]
%
% Outputs:
% nodeList - N(1)*N(2)x3 matrix with the index of each node in the 
%            first column, (x,y) position in the 2nd and 3rd
%            columns correspondigly
% elementList - Me*4 matrix. 
%               First column is the element index
%               2nd to 4th columns are the three vertices indexes
%               of the triangular element
% segmentList - Ms*3 matrix
%               First column is the segment index
%               2nd and 3rd columns are the starting point and
%               end point of the segment, correspondingly
%
% Meshing is done with the follwing form
% ________
% |     /|
% |    / |
% |   /  |
% |  /   |
% | /    |
% |/     |
% --------
%

dx = maxD(2);
dy = maxD(1);

% Generate the thee bottom meshes
obsLeft = obsCenter(2) - 0.5*obsDims(2);
obsRight = obsCenter(2) + 0.5*obsDims(2);
obsTop = obsCenter(1) + 0.5*obsDims(1);

% Bottom meshes
% bottom left mesh
blX = 0:dx:obsLeft;
% Obstacle mesh
bX = obsLeft:dx:obsRight;
bY = 0:dy:obsTop;
% bottom right mesh
brX = fliplr(rectDims(2):-dx:obsRight);

% Top mesh
tY = fliplr(rectDims(1):-dy:obsTop);

% Trunctuate edges to avoid double meshing
if(blX(end) == obsLeft) 
    blX = blX(1:end-1);
end
if(brX(1) == obsRight)
    brX = brX(2:end);
end
% Add the right edge if necessary
if(bX(end) ~= obsRight)
    bX = [bX obsRight];
end

if(bY(end) ~= obsTop)
    bY = [bY obsTop];
end
if(tY(1) == obsTop)
    tY = tY(2:end);
end

% Left to right mesh
lrX = [blX bX brX];
% Bottom to top mesh
btY = [bY tY];

[X,Y] = meshgrid(lrX,btY);

N = [numel(btY) numel(lrX)];

idxs = 1:N(1)*N(2);

nodeList = [idxs' X(:) Y(:)];

idxs = reshape(idxs,size(X));

% Map vertices of each squre
bottomLeft = idxs(1:end-1,1:end-1);
bottomRight = idxs(1:end-1,2:end);
topLeft = idxs(2:end,1:end-1);
topRight = idxs(2:end,2:end);

% Upper triangle list
uTriIdxs = 1:numel(topLeft);
elementListUpper = [uTriIdxs' topLeft(:) topRight(:) bottomLeft(:)];
lTriIdxs = numel(topLeft)+1:2*numel(topLeft);
elementListLower = [lTriIdxs' topRight(:) bottomRight(:) bottomLeft(:)];

elementList = [elementListUpper ; elementListLower];

% Now map surface segments
leftSide = idxs(1:end,1);
topSide = idxs(end,2:end);
rightSide = idxs(end-1:-1:1,end);
bottomSide = idxs(1,end-1:-1:2);

completeSideList = [leftSide(:) ; topSide(:) ; rightSide(:) ; bottomSide(:)];
numSegments = numel(completeSideList);
segmentList = [ (1:numSegments)' ...
                completeSideList ...
                [completeSideList(2:end) ; completeSideList(1)]];


