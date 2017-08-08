clc;
clear;
close all;

debug = 0;

e0 = (1e-9)*(36*pi)^-1;
m0 = 4*pi*1e-7;
c0 = 3e8;

% Amplitude of magnetic field 
H0 = 1;

a = 3.5e-2;
c = 5e-2;
f0 = 3e9;
wl = c0/f0; % Wave length
k0 = 2*pi*f0/c0;

hs = linspace(0.1*wl,0.35*wl,100);
R = zeros(size(hs));
T = zeros(size(hs));

myIdx = 1;
Zs = 10*(1+1i);

for h = hs;

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



%% Calculate K matrix and b vector

qs = zeros([size(segmentList,1) 1]);
qs(2:N(1)-1) = -2*1i*k0*H0;

gammas = zeros([size(segmentList,1) 1]);
% gammas(2:N(1)-1) = -1i*k0;
% gammas(N(1)+N(2):2*N(1)+N(2)) = -1i*k0;
gammas(2:N(1)-1) = -1i*k0;
gammas(N(1):N(1)+N(2)-1) = -1i*e0*Zs;
gammas(N(1)+N(2):2*N(1)+N(2)) = -1i*k0;
gammas(2*N(1)+N(2):end) = -1i*e0*Zs;

[K,b] = construct_Kmat_bvect(nodeList,elementList,segmentList,ax,ay,beta,gammas,qs);

%% Calculate Hz

Hz = K\b;

% Get indices of right\left edge
idxListLeft = find(nodeList(:,2) == 0);
idxListRight = find(nodeList(:,2) == rectSize(2));


HzRight = mean(Hz(idxListRight));
HzLeft = mean(Hz(idxListLeft));


R(myIdx) = (HzLeft-H0)/H0;

T(myIdx) = HzRight/(H0*exp(-1i*k0*rectSize(2)));

myIdx = myIdx + 1;

end

% figure;
% subplot(1,2,1);
% plot(hs/wl,abs(R).^2),xlabel('h/\lambda'),ylabel('|R|'),title('\epsilon_r = 4');
% subplot(1,2,2);
% plot(hs/wl,abs(T).^2),xlabel('h/\lambda'),ylabel('|T|'),title('\epsilon_r = 4');

figure;
% plot(hs/wl,abs(R).^2 + abs(T).^2)
plot(hs/wl,abs(R),'-b'),xlabel('h/\lambda'),ylabel('|R| or |T|'),title('\epsilon_r = 4');
hold on;
plot(hs/wl,abs(T),'-r');
hold off;
legend('|R|','|T|');

figure;
plot(hs/wl,abs(R).^2 + abs(T).^2),axis([min(hs/wl) max(hs/wl) 0 1.2]);
