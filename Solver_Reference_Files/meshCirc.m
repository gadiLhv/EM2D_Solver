function [nodeList,elementList,segmentList] = meshCirc(Rin,Rout,Nphi,Nr)

% Inputs
% R - Radius of circle
% Nphi - Number of divisions around the circle
% Nr - number of divisions along the line

rads = linspace(Rin,Rout,Nr+1);
phi = linspace(0,2*pi,Nphi+1);
phi = phi(1:end-1);

X = cos(phi).' * rads;
Y = sin(phi).' * rads;

nodeList = [(1:numel(X))' X(:) Y(:)];

idxSq = reshape(1:numel(X),size(X));

Ibl = idxSq(:,1:end-1);
Ibr = [idxSq(2:end,1:end-1) ; idxSq(1,1:end-1)];
Itl = idxSq(:,2:end);
Itr = [idxSq(2:end,2:end) ; idxSq(1,2:end)];

elementList = [(1:numel(Ibl)*2)'  [ ...
                Ibl(:) Ibr(:) Itl(:) ; ...
                Itl(:) Ibr(:) Itr(:) ]];
      
% Bottom boundary index list
if(Rin ~= 0)
    BBlistLeft = [];
    BBlistRight = [];
else
    BBlistLeft = idxSq(:,1);
    BBlistRight = [idxSq(2:end,1) ; idxSq(1,1)];
end

% Top boundary
TBlistLeft = idxSq(:,end);
TBlistRight = [idxSq(2:end,end) ; idxSq(1,end)];

segmentList = [(1:numel(TBlistLeft))' [BBlistLeft ; TBlistLeft] [BBlistRight ; TBlistRight]];