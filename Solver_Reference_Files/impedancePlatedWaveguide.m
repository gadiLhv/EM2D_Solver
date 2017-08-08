clc;
clear;
close all;

debug = 0;

e0 = (1e-9)*(36*pi)^-1;
m0 = 4*pi*1e-7;
c0 = 3e8;

a = 3.5e-2;
h = 1.85e-2;
c = 5e-2;
f0 = 3e9;
wl = c0/f0; % Wave length
k0 = 2*pi*f0/c0;

Zs = 10*(1+1i);

% Amplitude of magnetic field 
H0 = 1;

% Center of the obstacle in X,Y terms
oCenter = [h/2 1.5*c];
% Dimensions of obstacle
oDims = [h c];
% Relative permitivity of obstacle
er = 4;

%% Mesh the rectangle
rectSize = [a 3*c];

% [nodeList,elementList,segmentList] = meshRect(rectSize,N);
[nodeList,elementList,segmentList,N] = meshRect_wHobstacle(rectSize,[wl/40 wl/20],oDims,oCenter);

if(debug)
    figure;
    cVertIdxs = elementList(1,2:end);
    cVertCoors = nodeList(cVertIdxs,2:3);
    cColor = rand([1 3]);
    fill(cVertCoors(:,1),cVertCoors(:,2),cColor);
    axis([0 rectSize(2) 0 rectSize(1)]);
    xlabel('x'),ylabel('y');
    hold on;
    for i = 2:size(elementList,1)
        cVertIdxs = elementList(i,2:end);
        cVertCoors = nodeList(cVertIdxs,2:3);
        cColor = rand([1 3]);
        fill(cVertCoors(:,1),cVertCoors(:,2),cColor);
    end
    hold off;
    clear cVertCoors cVertIdxs i cColor
    
end

%% Assign specific alpha values for elements

% Calculate center of mass for each element
v = [elementList(:,2) elementList(:,3) elementList(:,4)];

CM = (nodeList(v(:,1),2:3) + nodeList(v(:,2),2:3) + nodeList(v(:,3),2:3))./3;

binArr =    (CM(:,1) >= (oCenter(2)-oDims(2)/2)) & ...
            (CM(:,2) >= (oCenter(1)-oDims(1)/2)) & ...
            (CM(:,1) <= (oCenter(2)+oDims(2)/2)) & ...
            (CM(:,2) <= (oCenter(1)+oDims(1)/2));

% Alpha arrays
ax = ones(size(binArr));
ay = ax;
ax(binArr) = 1/er;
ay(binArr) = 1/er;

beta = -k0.^2;

if(debug)
    figure;
    plot(0,0);
    axis([0 rectSize(2) 0 rectSize(1)]);
    xlabel('x'),ylabel('y');
    hold on;
    
    validElements = elementList(find(binArr),:);
    
    for i = 1:size(validElements,1)
        cVertIdxs = validElements(i,2:end);
        cVertCoors = nodeList(cVertIdxs,2:3);
        fill(cVertCoors(:,1),cVertCoors(:,2),'b');
        
        cCM = sum(cVertCoors)./3;
        plot(cCM(1),cCM(2),'.r');
        
    end
    
    plot([  oCenter(1)-oDims(1)/2 ...
            oCenter(1)-oDims(1)/2 ...
            oCenter(1)+oDims(1)/2 ...
            oCenter(1)+oDims(1)/2], ...
         [  oCenter(2)-oDims(2)/2 ...
            oCenter(2)+oDims(2)/2 ...
            oCenter(2)+oDims(2)/2 ...
            oCenter(2)-oDims(2)/2], ...
            '--g');
    hold off;
    
    clear cVertIdxs validElements i cVertCoors cCM
end

%% Calculate K matrix and b vector

qs = zeros([size(segmentList,1) 1]);
qs(2:N(1)-1) = -2*1i*k0*H0;

gammas = zeros([size(segmentList,1) 1]);
gammas(2:N(1)-1) = -1i*k0;
gammas(N(1):N(1)+N(2)-1) = -1i*e0*Zs;
gammas(N(1)+N(2):2*N(1)+N(2)) = -1i*k0;
gammas(2*N(1)+N(2):end) = -1i*e0*Zs;
[K,b] = construct_Kmat_bvect(nodeList,elementList,segmentList,ax,ay,beta,gammas,qs);

%% Calculate Hz

Hz = K\b;


%% Thicken the mesh and draw
nMesh = [200 100];
[xMesh,yMesh,HzI] = thickenMeshByElement(nodeList,elementList,Hz,rectSize,nMesh);

% Calculate electric fields
dx = xMesh(1,2)-xMesh(1,1);
dy = yMesh(2,1)-yMesh(1,1);
ExI = zeros(size(HzI));
EyI = ExI;

ExI(2:end-1,:) = 1i*2*pi*f0*e0*(HzI(3:end,:) - HzI(1:end-2,:))/(2*dy);
EyI(:,2:end-1) = -1i*2*pi*f0*e0*(HzI(:,3:end) - HzI(:,1:end-2))/(2*dx);

figure;
subplot(3,1,1);
hdl = pcolor(xMesh,yMesh,real(HzI));
set(hdl,'linestyle','none'),xlabel('x'),ylabel('y'),title('H_z');
hold on
plot([  oCenter(1)-oDims(1)/2 oCenter(1)-oDims(1)/2 oCenter(1)+oDims(1)/2 oCenter(1)+oDims(1)/2], ...
     [  oCenter(2)-oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)-oDims(2)/2], ...
        '-r','linewidth',2);
hold off;

subplot(3,1,2);
hdl = pcolor(xMesh,yMesh,real(ExI));
set(hdl,'linestyle','none'),xlabel('x'),ylabel('y'),title('E_x');
hold on
plot([  oCenter(1)-oDims(1)/2 oCenter(1)-oDims(1)/2 oCenter(1)+oDims(1)/2 oCenter(1)+oDims(1)/2], ...
     [  oCenter(2)-oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)-oDims(2)/2], ...
        '-r','linewidth',2);
hold off;

subplot(3,1,3);
hdl = pcolor(xMesh,yMesh,real(EyI));
set(hdl,'linestyle','none'),xlabel('x'),ylabel('y'),title('E_y');
hold on
plot([  oCenter(1)-oDims(1)/2 oCenter(1)-oDims(1)/2 oCenter(1)+oDims(1)/2 oCenter(1)+oDims(1)/2], ...
     [  oCenter(2)-oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)-oDims(2)/2], ...
        '-r','linewidth',2);
hold off;

