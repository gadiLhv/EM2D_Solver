clc;
clear;
close all;

debug = 1;

e0 = (1e-9)*(36*pi)^-1;
m0 = 4*pi*1e-7;
c0 = 3e8;

a = 3.5e-2;
h = 1.85e-2;
c = 5e-2;
f0 = 3e9;
wl = c0/f0; % Wave length
k0 = 2*pi*f0/c0;

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
    patch('faces',elementList(:,2:end),'vertices',nodeList(:,2:3), ...
          'facecolor',[0.3,0.8,0.3], ...
          'edgecolor',[0,0,0]) ;

    axis([0 rectSize(2) 0 rectSize(1)]);
    xlabel('x'),ylabel('y');
    
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
    
    patch('faces',elementList(:,2:end),'vertices',nodeList(:,2:3), ...
      'facecolor',[0.3,0.8,0.3], ...
      'edgecolor',[0,0,0]) ;
    
    axis([0 rectSize(2) 0 rectSize(1)]);
    xlabel('x'),ylabel('y');
    
    
    hold on;
    validElements = elementList(find(binArr),2:end);
    
    patch('faces',validElements,'vertices',nodeList(:,2:3), ...
      'facecolor',[0.3,0.3,0.8], ...
      'edgecolor',[0,0,0]) ;
    
    plot([  oCenter(1)-oDims(1)/2 ...
            oCenter(1)-oDims(1)/2 ...
            oCenter(1)+oDims(1)/2 ...
            oCenter(1)+oDims(1)/2], ...
         [  oCenter(2)-oDims(2)/2 ...
            oCenter(2)+oDims(2)/2 ...
            oCenter(2)+oDims(2)/2 ...
            oCenter(2)-oDims(2)/2], ...
            '--g');
    
end

%% Calculate K matrix and b vector

qs = zeros([size(segmentList,1) 1]);
qs(2:N(1)-1) = -2*1i*k0*H0;

gammas = zeros([size(segmentList,1) 1]);
gammas(2:N(1)-1) = -1i*k0;
gammas(N(1)+N(2):2*N(1)+N(2)) = -1i*k0;

[K,b] = construct_Kmat_bvect(nodeList,elementList,segmentList,ax,ay,beta,gammas,qs);

%% Calculate Hz

Hz = K\b;


%% Thicken the mesh and draw
nMesh = [200 100];
[xMesh,yMesh,HzI] = thickenMeshByElement(nodeList,elementList,Hz,[0 rectSize(2) 0 rectSize(1)],nMesh);

% Calculate electric fields
dx = xMesh(1,2)-xMesh(1,1);
dy = yMesh(2,1)-yMesh(1,1);
ExI = zeros(size(HzI));
EyI = ExI;

ExI(2:end-1,:) = 1i*2*pi*f0*e0*(HzI(3:end,:) - HzI(1:end-2,:))/(2*dy);
EyI(:,2:end-1) = -1i*2*pi*f0*e0*(HzI(:,3:end) - HzI(:,1:end-2))/(2*dx);

figure;
subplot(3,1,1);
% hdl = pcolor(xMesh,yMesh,real(HzI));
hdl = contourf(xMesh,yMesh,real(HzI));
% set(hdl,'linestyle','none');
xlabel('x'),ylabel('y'),title('H_z');
hold on
plot([  oCenter(1)-oDims(1)/2 oCenter(1)-oDims(1)/2 oCenter(1)+oDims(1)/2 oCenter(1)+oDims(1)/2], ...
     [  oCenter(2)-oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)-oDims(2)/2], ...
        '-r','linewidth',3);
hold off;

subplot(3,1,2);
% hdl = pcolor(xMesh,yMesh,real(ExI));
hdl = contourf(xMesh,yMesh,real(ExI));
% set(hdl,'linestyle','none');
xlabel('x'),ylabel('y'),title('E_x');
hold on
plot([  oCenter(1)-oDims(1)/2 oCenter(1)-oDims(1)/2 oCenter(1)+oDims(1)/2 oCenter(1)+oDims(1)/2], ...
     [  oCenter(2)-oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)-oDims(2)/2], ...
        '-r','linewidth',3);
hold off;

subplot(3,1,3);
% hdl = pcolor(xMesh,yMesh,real(EyI));
hdl = contourf(xMesh,yMesh,real(EyI));
% set(hdl,'linestyle','none');
xlabel('x'),ylabel('y'),title('E_y');
hold on
plot([  oCenter(1)-oDims(1)/2 oCenter(1)-oDims(1)/2 oCenter(1)+oDims(1)/2 oCenter(1)+oDims(1)/2], ...
     [  oCenter(2)-oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)+oDims(2)/2 oCenter(2)-oDims(2)/2], ...
        '-r','linewidth',3);
hold off;

