function [nodeList,elementList,segmentList] = meshCircInAndOut(Rhollow,Rd,Rout,Nphi,Nr)

% Inputs
% R - Radius of circle
% Nphi - Number of divisions around the circle
% Nr - number of divisions along the line

% Radii of the hollow cylinder
rads1 = linspace(0,Rhollow,Nr+1);
rads1 = rads1(2:end);
% Radii of the dielectric cylinder
rads2 = linspace(Rhollow,Rd,Nr+1);
rads2 = rads2(2:end);
% Radii of the outer area
rads3 = linspace(Rd,Rout,Nr+1);
rads3 = rads3(2:end);

phi = linspace(0,2*pi,Nphi+1);
phi = phi(1:end-1);

X1 = cos(phi).' * rads1;
Y1 = sin(phi).' * rads1;

X2 = cos(phi).' * rads2;
Y2 = sin(phi).' * rads2;

X3 = cos(phi).' * rads3;
Y3 = sin(phi).' * rads3;

X = [X1 X2 X3];
Y = [Y1 Y2 Y3];

nodeList = [(1:numel(X))' X(:) Y(:)];

idxSq = reshape(1:numel(X),size(X));

Ibl = idxSq(:,1:end-1);
Ibr = [idxSq(2:end,1:end-1) ; idxSq(1,1:end-1)];
Itl = idxSq(:,2:end);
Itr = [idxSq(2:end,2:end) ; idxSq(1,2:end)];

elementList = [(1:numel(Ibl)*2)'  [ ...
                Ibl(:) Ibr(:) Itl(:) ; ...
                Itl(:) Ibr(:) Itr(:) ]];   

% Add the center as the last node
nodeList = [nodeList ; [numel(X)+1 0 0]];

% Construct the additional triangles
leftIdxs = idxSq(:,1);
rightIdxs = [idxSq(2:end,1) ; idxSq(1,1)];
eAdd = [(numel(Ibl)*2+1:numel(Ibl)*2+Nphi)' ...
        ones([Nphi 1])*size(nodeList,1) ...
        leftIdxs(:) ...
        rightIdxs(:) ] ;
elementList = [ elementList ; eAdd ];      

% Bottom boundary index list
% if(Rin == 0)
BBlistLeft = [];
BBlistRight = [];
% else
%     BBlistLeft = idxSq(:,1);
%     BBlistRight = [idxSq(2:end,1) ; idxSq(1,1)];
% end

% Top boundary
TBlistLeft = idxSq(:,end);
TBlistRight = [idxSq(2:end,end) ; idxSq(1,end)];

segListLeft = [BBlistLeft ; TBlistLeft];
segListRight = [BBlistRight ; TBlistRight];

segmentList = [(1:numel(segListLeft))' segListLeft segListRight];