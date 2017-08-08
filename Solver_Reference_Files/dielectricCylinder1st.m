clc
clear
close all

debug = 0;

e0 = (1e-9)*(36*pi)^-1;
m0 = 4*pi*1e-7;
c0 = 3e8;

f0 = 3e9;
w0 = 2*pi*f0;
wl = c0/f0;
k0 = 2*pi/wl;

H0 = 1;

er = 4;

a = 0.25*wl;
b = 0.3*wl;
c = 0.4*wl;

Rin = a;
Rout = b;
Rabs = c;

Nr = 3;
Nphi = 44;

% Impact angle
phi0 = pi;

%%

[nodeList,elementList,segmentList] = meshCircAndOut(Rin,Rout,Rabs,Nphi,Nr);

if(debug)
    figure;
    cVertIdxs = elementList(1,2:end);
    cVertCoors = nodeList(cVertIdxs,2:3);
    cColor = rand([1 3]);
    fill(cVertCoors(:,1),cVertCoors(:,2),cColor,'edgecolor','none');
    axis([-c c -c c]);
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
end
%% Assign specific alpha values for elements

v = [elementList(:,2) elementList(:,3) elementList(:,4)];

% Calculate center of mass for each element
CM = (nodeList(v(:,1),2:3) + nodeList(v(:,2),2:3) + nodeList(v(:,3),2:3))./3;

binArr = (sqrt(CM(:,1).^2 + CM(:,2).^2) >= Rin) & (sqrt(CM(:,1).^2 + CM(:,2).^2) <= Rout);

% Alpha arrays
ax = ones(size(binArr));
ay = ax;
ax(binArr) = 1/er;
ay(binArr) = 1/er;

beta = -k0.^2;

if(debug)
    figure;
    plot(0,0);
    axis([-c c -c c]);
    xlabel('x'),ylabel('y');
    hold on;
    
    validElements = elementList(find(binArr),:);
    
    for i = 1:size(validElements,1)
        cVertIdxs = validElements(i,2:end);
        cVertCoors = nodeList(cVertIdxs,2:3);
        fill(cVertCoors(:,1),cVertCoors(:,2),'b','edgecolor','none');
        
        cCM = sum(cVertCoors)./3;
        plot(cCM(1),cCM(2),'.r');
        
    end
    
    phiTemp = linspace(0,2*pi,100);
    
    plot(Rin*cos(phiTemp),Rin*sin(phiTemp),'--g');
    plot(Rout*cos(phiTemp),Rout*sin(phiTemp),'--g');
    hold off;
    
    clear cVertIdxs validElements i cVertCoors cCM phiTemp
end

%% Calulcate approprate gamma and q

segmentNodeList = segmentList(:,2:3);

% Get coordinates of segment nodes
sNodeCoorsLeft = nodeList(segmentNodeList(:,1),2:3);
sNodeCoorsRight = nodeList(segmentNodeList(:,2),2:3);
% Mean between left and right
sNodeCoors = (sNodeCoorsLeft + sNodeCoorsRight)./2;
% The angle is allready the mean value
sNodeAngles = atan2(sNodeCoors(:,2),sNodeCoors(:,1));

% Get coordinates of outermost nodes for source infusion
% Do this by checking if the positions is outside the middle ring
outerBin = (sqrt(sNodeCoorsLeft(:,1).^2 + sNodeCoorsLeft(:,2).^2) >= Rout) & ...
           (sqrt(sNodeCoorsRight(:,1).^2 + sNodeCoorsRight(:,2).^2) >= Rout);

q = zeros([size(segmentList,1) 1]);
q(outerBin) = exp(1i*k0*Rabs*cos(sNodeAngles(outerBin)-phi0)).* ...
                (1i*k0*cos(sNodeAngles(outerBin)-phi0) + 1i*k0 + 1./(2*Rabs));


% Now gammas
gamma = zeros([size(segmentList,1) 1]);
gamma(outerBin) = 1i*k0 + 1./(2*Rabs);

%% Solve!!!
[K,b] = construct_Kmat_bvect(nodeList,elementList,segmentList,ax,ay,beta,gamma,q);
Hz = K\b;

%% Thicken the mesh and draw
nMesh = [400 400];
[xMesh,yMesh,HzI] = thickenMeshByElement(nodeList,elementList,Hz,[-c c -c c],nMesh);

% Add the incident field
HzI = HzI + exp(1i*k0*(xMesh*cos(phi0)+yMesh*sin(phi0)));

% Calculate electric fields
dx = xMesh(1,2)-xMesh(1,1);
dy = yMesh(2,1)-yMesh(1,1);
ExI = zeros(size(HzI));
EyI = ExI;

ExI(2:end-1,:) = 1i*2*pi*f0*e0*(HzI(3:end,:) - HzI(1:end-2,:))/(2*dy);
EyI(:,2:end-1) = -1i*2*pi*f0*e0*(HzI(:,3:end) - HzI(:,1:end-2))/(2*dx);

figure;
subplot(1,3,1);
hdl = pcolor(xMesh,yMesh,abs(HzI));
% hdl = contourf(xMesh,yMesh,real(HzI));
set(hdl,'linestyle','none');
xlabel('x'),ylabel('y'),title('H_z');
hold on
phiTemp = linspace(0,2*pi,100);
    
plot(Rin*cos(phiTemp),Rin*sin(phiTemp),'--r');
plot(Rout*cos(phiTemp),Rout*sin(phiTemp),'--r');
hold off;

subplot(1,3,2);
hdl = pcolor(xMesh,yMesh,abs(ExI));
% hdl = contourf(xMesh,yMesh,real(ExI));
set(hdl,'linestyle','none');
xlabel('x'),ylabel('y'),title('E_x');
hold on
phiTemp = linspace(0,2*pi,100);
    
plot(Rin*cos(phiTemp),Rin*sin(phiTemp),'--r');
plot(Rout*cos(phiTemp),Rout*sin(phiTemp),'--r');
hold off;

subplot(1,3,3);
hdl = pcolor(xMesh,yMesh,abs(EyI));
% hdl = contourf(xMesh,yMesh,real(EyI));
set(hdl,'linestyle','none');
xlabel('x'),ylabel('y'),title('E_y');
hold on
phiTemp = linspace(0,2*pi,100);
    
plot(Rin*cos(phiTemp),Rin*sin(phiTemp),'--r');
plot(Rout*cos(phiTemp),Rout*sin(phiTemp),'--r');
hold off;

%% Draw the inner current

phiTag = linspace(0,pi,100);

Xc = (Rin + 0.00001)*cos(phiTag);
Yc = (Rin + 0.00001)*sin(phiTag);

HzsC = fillMeshByElement(nodeList,elementList,Hz,Xc,Yc);

if(sum(isnan(HzsC(:))) ~= 0) error('Assertion error'); end;

HzC = HzsC + exp(1i*k0*(Xc.*cos(phi0)+Yc*sin(phi0)));
Hzi = exp(1i*k0*(Xc.*cos(phi0)+Yc*sin(phi0)));

figure;
plot(phiTag,abs(HzC./Hzi))