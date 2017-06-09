% test_mesh.m

clc
clear
close all

mesherPath = './Mesh2d_v23';

%%


addpath(mesherPath);

%% Set options for solver

% Set printout/no printout
outputMesh = false(1);

% Global mesh properties

% 1. Max volume step for background
maxBackgrondStep = 0.4;

% 2. Max step for amaterials
maxMaterialStep = 0.1;


% Local mesh properties
% 

%% Problem dimensions
% Air cube dimensions (background)
boundingBox = [ 0 0 
                4 0
                4 4
                0 4];
                
boundingBoxConnectivity = [ 1 2
                            2 3
                            3 4
                            4 1];
                            

% Mesh properties for background

                            
                            
%% Obstacle definitions

% Rectangular scatterer
rectangularScatterer = [0.5 0.5
                        0.5 2.5
                        2.5   2.5
                        2.5   0.5];
rectangularScattererConnectivity = [1 2
                                    2 3
                                    3 4
                                    4 1];



                                    

% Elliptic scatterer
Npts = 40;
ellipticRad = 0.65;
ellipticSkew = [1 0.55]; 
ellipticRotate = 30*180/pi;
ellipseCenter = [2.2 2.2];
thetaPoints = linspace(0,2*pi,Npts+1);
thetaPoints = thetaPoints(1:end-1);
ellipticScatterer = ellipticRad*[cos(thetaPoints.')*ellipticSkew(1) sin(thetaPoints.')*ellipticSkew(2)];
ellipticScatterer  = ([cos(ellipticRotate) -sin(ellipticRotate) ; sin(ellipticRotate) cos(ellipticRotate)]*(ellipticScatterer.')).';
ellipticScatterer = bsxfun(@plus,ellipticScatterer,ellipseCenter);
ellipticScattererConnectivity = 1:size(ellipticScatterer,1);
ellipticScattererConnectivity = [ellipticScattererConnectivity.' [ellipticScattererConnectivity(2:end).' ; 1]];


%% Test mesher

% Cascade nodes
%nodes = [boundingBox          ; ...
%         rectangularScatterer ; ...
%         ellipticScatterer];
nodes = [boundingBox          ; ...
         rectangularScatterer];



% Define edges
%connectivity = [boundingBoxConnectivity ; ... 
%                rectangularScattererConnectivity+size(boundingBoxConnectivity,1) ; ...
%                ellipticScattererConnectivity+size(boundingBoxConnectivity,1)+size(rectangularScattererConnectivity,1)];
connectivity = [boundingBoxConnectivity ; ... 
                rectangularScattererConnectivity+size(boundingBoxConnectivity,1)];


% Force meshing for each shape separately by assigning edges to faces
%faces = cell([3 1]);
%faces{1} = 1:size(connectivity,1);
%faces{2} = (1:size(rectangularScattererConnectivity,1))+size(boundingBoxConnectivity,1);
%faces{3} = (1:size(ellipticScattererConnectivity,1))+size(boundingBoxConnectivity,1)+size(rectangularScattererConnectivity,1);
faces = cell([2 1]);
faces{1} = 1:size(connectivity,1);
faces{2} = (1:size(rectangularScattererConnectivity,1))+size(boundingBoxConnectivity,1);


                
% Define mesh properties 
%hdata.edgeh = [faces{1}.' maxBackgrondStep*ones(size(faces{1})).' ; ...
%               faces{2}.' maxMaterialStep*ones(size(faces{2})).' ; ...
%               faces{3}.' maxMaterialStep*ones(size(faces{3})).'];
hdata.edgeh = [faces{1}.' maxBackgrondStep*ones(size(faces{1})).' ; ...
               faces{2}.' maxMaterialStep*ones(size(faces{2})).' ];



                
[p,t] = meshfaces(nodes,connectivity,faces,hdata);

               
%rmpath(mesherPath);