clc
clear
close all

f0 = 3e9;
c0 = 3e8;

wl = c0/f0;

% Rin = 0.25*wl;
% Rout = 0.3*wl;
Rin = 0.5;
Rout = 1.2;
Rabs = 1.8;

Nr = 3;
Nphi = 44;

% [nodeList,elementList,segmentList] = meshCirc(Rin,Rout,Nphi,Nr);
% [nodeList,elementList,segmentList] = meshCircAndOut(Rin,Rout,Rabs,Nphi,Nr);
[nodeList,elementList,segmentList] = meshCircInAndOut(Rin,Rout,Rabs,Nphi,Nr);



figure;
cVertIdxs = elementList(1,2:end);
cVertCoors = nodeList(cVertIdxs,2:3);
cColor = rand([1 3]);
h = fill(cVertCoors(:,1),cVertCoors(:,2),cColor,'edgecolor','none');
axis([-2 2 -2 2]);
xlabel('x'),ylabel('y');
hold on;
for i = 2:size(elementList,1)
    cVertIdxs = elementList(i,2:end);
    cVertCoors = nodeList(cVertIdxs,2:3);
    cColor = rand([1 3]);
    fill(cVertCoors(:,1),cVertCoors(:,2),cColor,'edgecolor','none');
end
hold off;
clear cVertCoors cVertIdxs i cColor
