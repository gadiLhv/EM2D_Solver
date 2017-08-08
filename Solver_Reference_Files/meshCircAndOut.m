function [nodeList,elementList,segmentList] = meshCircAndOut(Rin,Rout,Rabs,Nphi,Nr)

% Inputs
% R - Radius of circle
% Nphi - Number of divisions around the circle
% Nr - number of divisions along the line

rads1 = linspace(Rin,Rout,Nr+1);
rads2 = linspace(Rout,Rabs,Nr+1);
rads2 = rads2(2:end);
% Exclide middle point if there is no hollow cylinder
if(Rin == 0)
    rads1 = rads1(2:end);
end
phi = linspace(0,2*pi,Nphi+1);
phi = phi(1:end-1);

X1 = cos(phi).' * rads1;
Y1 = sin(phi).' * rads1;

X2 = cos(phi).' * rads2;
Y2 = sin(phi).' * rads2;

X = [X1 X2];
Y = [Y1 Y2];

nodeList = [(1:numel(X))' X(:) Y(:)];

idxSq = reshape(1:numel(X),size(X));

Ibl = idxSq(:,1:end-1);
Ibr = [idxSq(2:end,1:end-1) ; idxSq(1,1:end-1)];
Itl = idxSq(:,2:end);
Itr = [idxSq(2:end,2:end) ; idxSq(1,2:end)];

elementList = [(1:numel(Ibl)*2)'  [ ...
                Ibl(:) Ibr(:) Itl(:) ; ...
                Itl(:) Ibr(:) Itr(:) ]];   

% If there is no hollow, construct a special addition
if(Rin == 0)
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
end       
% Bottom boundary index list
if(Rin == 0)
    BBlistLeft = [];
    BBlistRight = [];
else
    BBlistLeft = idxSq(:,1);
    BBlistRight = [idxSq(2:end,1) ; idxSq(1,1)];
end

% Top boundary
TBlistLeft = idxSq(:,end);
TBlistRight = [idxSq(2:end,end) ; idxSq(1,end)];

segListLeft = [BBlistLeft ; TBlistLeft];
segListRight = [BBlistRight ; TBlistRight];

segmentList = [(1:numel(segListLeft))' segListLeft segListRight];