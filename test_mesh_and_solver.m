% test_mesh.m

clc
clear
close all

addpath("./FEM2Dsolver_V01");
addpath("./Mesh2d_v23");

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
                3 0
                3 3
                0 3];
                
boundingBoxConnectivity = [ 1 2
                            2 3
                            3 4
                            4 1];
                            

% Mesh properties for background

                            
                            
%% Obstacle definitions

% Rectangular scatterer
rectangularScatterer = [0.5 0.5
                        0.5 1
                        1   1
                        1   0.5];
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
nodes = [boundingBox
         rectangularScatterer
         ellipticScatterer];

% Define edges
connectivity = [boundingBoxConnectivity
                rectangularScattererConnectivity+size(boundingBoxConnectivity,1)
                ellipticScattererConnectivity+size(boundingBoxConnectivity,1)+size(rectangularScattererConnectivity,1)
                ];

% Force meshing for each shape separately by assigning edges to faces
faces = cell([3 1]);
faces{1} = 1:size(connectivity,1);
faces{2} = (1:size(rectangularScattererConnectivity,1))+size(boundingBoxConnectivity,1);
faces{3} = (1:size(ellipticScattererConnectivity,1))+size(boundingBoxConnectivity,1)+size(rectangularScattererConnectivity,1);
                
% Define mesh properties 
hdata.edgeh = [faces{1}.' maxBackgrondStep*ones(size(faces{1})).'
               faces{2}.' maxMaterialStep*ones(size(faces{2})).'
               faces{3}.' maxMaterialStep*ones(size(faces{3})).'];

  
options.output = outputMesh;
                
% Mesh the mofo!
[points,tetraNodes,faceAssignment,edgeAssignments] = meshfaces(nodes,connectivity,faces,hdata,options);

figure;
plot(0,0);
axis([0 3 0 3]);
hold on;
shapecolors = 'rgb';
%% Prepare materials for solver
%for faceIdx = 1:numel(faces)
%  % Get corresponding tetrahedrons for each face
%  tetraIdxs = find(faceAssignment == faceIdx);
%  cNodes = tetraNodes(tetraIdxs,:);
%  cColor = shapecolors(mod(faceIdx-1,numel(shapecolors))+1);
%  for nodeIdx = 1:size(cNodes,1)
%    % Get this triangle's points
%    cPts = points(cNodes(nodeIdx,:).',:);
%    plot([cPts(:,1) ; cPts(1,1)],[cPts(:,2) ; cPts(1,2)],['-' cColor]);
%  end 
%end
%
for faceIdx = 1:numel(edgeAssignments)
  % Get point indexes corresponding to edges
  cEdges = edgeAssignments{faceIdx};
  % And the two limits of each edge
  cEdgeCoors1 = points(cEdges(:,1),:);
  cEdgeCoors2 = points(cEdges(:,2),:);
  
  % Choose color
  cColor = shapecolors(mod(faceIdx-1,numel(shapecolors))+1);
  % Plot!!
  for eIdx = 1:size(cEdges,1)
    plot([cEdgeCoors1(eIdx,1) cEdgeCoors2(eIdx,1)],...
         [cEdgeCoors1(eIdx,2) cEdgeCoors2(eIdx,2)],...
         ['-' cColor],'linewidth',4);
  end
 
end
hold off;
