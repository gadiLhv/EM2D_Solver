function [nodeList,elementList,segmentList] = meshCircSq(R,N)

% R = 1;
% N = 10;

% Inputs
% R - Radius of circle
% N - number of samples between 0 and pi/2

% N = N;

dphi = (pi/2)/N;
phi = (0:N-1)*dphi;

xs = R*cos(phi);
ys = R*sin(phi);

xs = [-xs 0 fliplr(xs)];
ys = [-fliplr(ys) 0 ys]; 


% xs = sort(unique(xs));
% ys = sort(unique(ys));

[xMesh,yMesh] = meshgrid(xs,ys);

nodeList = [(1:numel(xMesh))' xMesh(:) yMesh(:)];

nodeMesh = reshape((1:numel(xMesh))',[numel(ys) numel(xs)]);

% Get indexes of all four vertices of each square
botLeft = nodeMesh(1:end-1,1:end-1);
botRight = nodeMesh(1:end-1,2:end);
topRight = nodeMesh(2:end,2:end);
topLeft = nodeMesh(2:end,1:end-1);

% q1 and q3 use main diagonal, q2 and q4 use secondary diagonal
q1BotLeft = botLeft(end/2+1:end,end/2+1:end);
q2BotLeft = botLeft(end/2+1:end,1:end/2);
q3BotLeft = botLeft(1:end/2,1:end/2);
q4BotLeft = botLeft(1:end/2,end/2+1:end);

q1BotRight = botRight(end/2+1:end,end/2+1:end);
q2BotRight = botRight(end/2+1:end,1:end/2);
q3BotRight = botRight(1:end/2,1:end/2);
q4BotRight = botRight(1:end/2,end/2+1:end);

q1TopRight = topRight(end/2+1:end,end/2+1:end);
q2TopRight = topRight(end/2+1:end,1:end/2);
q3TopRight = topRight(1:end/2,1:end/2);
q4TopRight = topRight(1:end/2,end/2+1:end);

q1TopLeft = topLeft(end/2+1:end,end/2+1:end);
q2TopLeft = topLeft(end/2+1:end,1:end/2);
q3TopLeft = topLeft(1:end/2,1:end/2);
q4TopLeft = topLeft(1:end/2,end/2+1:end);

% For q1 and q3, create upper right and lower left triangles
ur = [  [q1TopLeft(:) q1TopRight(:) q1BotRight(:)] ; ...
        [q3TopLeft(:) q3TopRight(:) q3BotRight(:)]];
ll = [  [q1TopLeft(:) q1BotRight(:) q1BotLeft(:)] ; ...
        [q3TopLeft(:) q3BotRight(:) q3BotLeft(:)]];

% For q2 and q4, create upper left and lower right triangles
ul = [  [q2TopLeft(:) q2TopRight(:) q2BotLeft(:)] ; ...
        [q4TopLeft(:) q4TopRight(:) q4BotLeft(:)]];
lr = [  [q2TopRight(:) q2BotRight(:) q2BotLeft(:)] ; ...
        [q4TopRight(:) q4BotRight(:) q4BotLeft(:)]];

elementList = [(1:size(ul,1)*4)' [ur ; ll ; ul ; lr]];

% Now map surface segments
leftSide = nodeMesh(1:end,1);
topSide = nodeMesh(end,2:end);
rightSide = nodeMesh(end-1:-1:1,end);
bottomSide = nodeMesh(1,end-1:-1:2);

completeSideList = [leftSide(:) ; topSide(:) ; rightSide(:) ; bottomSide(:)];
numSegments = numel(completeSideList);
segmentList = [ (1:numSegments)' ...
                completeSideList ...
                [completeSideList(2:end) ; completeSideList(1)]];

