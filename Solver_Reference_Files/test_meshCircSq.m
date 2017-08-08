clc
clear
close all

R = 1.8;
N = 12;

[nodeList , elementList , segmentList] = meshCircSq(R,N);

%%
figure;
cVertIdxs = elementList(1,2:end);
cVertCoors = nodeList(cVertIdxs,2:3);
cColor = rand([1 3]);
fill(cVertCoors(:,1),cVertCoors(:,2),cColor);
axis([-2 2 -2 2]);
xlabel('x'),ylabel('y');
hold on;
for i = 2:size(elementList,1)
    cVertIdxs = elementList(i,2:end);
    cVertCoors = nodeList(cVertIdxs,2:3);
    cColor = rand([1 3]);
    fill(cVertCoors(:,1),cVertCoors(:,2),cColor);
end
phi = linspace(0,2*pi,100);
plot(R*cos(phi),R*sin(phi),'-.r','linewidth',3);
hold off;
clear cVertCoors cVertIdxs i cColor