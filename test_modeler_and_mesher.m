% test_modeler_and_mesher
clc
clear
close all

modelerPath = './Modeler_2D';
mesherPath = './Mesh2D_Latest';

%% Mesher parameters
% Set printout/no printout
outputMesh = true(1);

% Global mesh properties

% 1. Max volume step for background
maxBackgrondStep = 1.5;

% 2. Max step for amaterials
maxMaterialStep = 0.6;

%% Define background
bgBoundX = [-4 4];
bgBoundY = [-4 4];

%% Model this shit!

addpath(modelerPath);
addpath(mesherPath);

bg = mod2D_createRectangleStruct([bgBoundX(1) bgBoundY(1)],[bgBoundX(2) bgBoundY(2)]);

% Geometry 1
pol1 = mod2D_createRectangleStruct([-2 -2],[2 -1]);
pol2 = mod2D_createRectangleStruct([-2 -1],[2 0]);
pol3 = mod2D_createRectangleStruct([-2 0],[2 1]);
rect1 = mod2D_createRectangleStruct([-1 -1.5],[1 0.5]);
pol1 =  mod2D_booleanOperation(pol1,rect1,'subtract');
pol2 =  mod2D_booleanOperation(pol2,rect1,'subtract');
pol3 =  mod2D_booleanOperation(pol3,rect1,'subtract');

% Geometery 2
%rect1 = mod2D_createRectangleStruct([-3 -3],[3 3]);
%rect2 = mod2D_createRectangleStruct([-2.5 -2.5],[2.5 2.5]);
%pol1 = mod2D_booleanOperation(rect1,rect2,'subtract');
%rect1 = mod2D_createRectangleStruct([-2 -2],[2 2]);
%rect2 = mod2D_createRectangleStruct([-1 -1],[1 1]);
%pol2 = mod2D_booleanOperation(rect1,rect2,'subtract');
%pol1 = mod2D_booleanOperation(pol1,pol2,'add');
%rect1 = mod2D_createRectangleStruct([2 -0.5],[2.5 1.5]);
%pol1 = mod2D_booleanOperation(pol1,rect1,'add');
%pol2 = mod2D_createRectangleStruct([-2.5 1],[-2 2]); % Bad at handling duplicate nodes. Remove them while creating connectivty list
%pol3 = mod2D_createPolygonStruct([-0.7 0.4 1],[-0.6 -0.8 1]);


polList = {pol1,pol2,pol3};

% Add all polygons before subtracting from the background
if(numel(polList) < 2)
  allPols = polList{1};
elseif(numel(polList) >= 2)
  allPols = mod2D_booleanOperation(polList{1},polList{2},'add');
  if(numel(polList) > 2)
    for polIdx = 3:numel(polList)
      allPols = mod2D_booleanOperation(allPols,polList{polIdx},'add');
    end
  end
end

% Subtract all polygons from background
bg = mod2D_booleanOperation(bg,allPols,'subtract');

polList = [{bg}  polList];
[node,edge,face] = mod2D_polygonToFaceList(polList,1e-8);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial feature size estimation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hmax = maxBackgrondStep;

% Estimate Local Feature Size (LFS) for each node\part
[vlfs,tlfs,hlfs] = lfshfn2( node,...
                            edge,...
                            face);

                            
% Where the length come out larger than hmax, limit.
hlfs = min(hmax,hlfs) ;

% Super special indexing function for AABB tree queries
slfs = idxtri2(vlfs,tlfs);

% Something to do with prior triangular indexing
hfun = @trihfn2;

% Refine initial mesh
[ vert,...              % Vertices of mesh cells
  etri,...              % Constrained edge list (???)
  tria,...              % Triangle threesomes (attached to VERT)
  tnum] = ...           % Part assignments
  refine2(node, ...     % Full node list
          edge, ...     % Edge connectivity
          face, ...     % Edge->Part assignment
          [], ...       % Options field (RHO not used here. Test!)
          hfun, ...     % Cool meshing limiting function
          vlfs,...      % Arguments for this function
          tlfs,...
          slfs,...
          hlfs) ;


%
figure;
subplot(1,2,1);
patch('faces',tria(tnum==1,1:3),'vertices',vert, ...
    'facecolor',[1.,1.,1.], ...
    'edgecolor',[0,0,0]) ;
hold on; 
axis image off;
randColor = rand([1 3]);
patch('faces',tria(tnum==2,1:3),'vertices',vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
randColor = rand([1 3]);
patch('faces',tria(tnum==3,1:3),'vertices',vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
randColor = rand([1 3]);
patch('faces',tria(tnum==4,1:3),'vertices',vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
title(['MESH.: KIND=DELFRONT, |TRIA|=', ...
    num2str(size(tria,1))]) ;
hold off;

% Smooth mesh
[vert,etri,tria,tnum] = smooth2(vert,etri,tria,tnum) ;

subplot(1,2,2);
patch('faces',tria(tnum==1,1:3),'vertices',vert, ...
    'facecolor',[1.,1.,1.], ...
    'edgecolor',[0,0,0]) ;
hold on; 
axis image off;
randColor = rand([1 3]);
patch('faces',tria(tnum==2,1:3),'vertices',vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
randColor = rand([1 3]);
patch('faces',tria(tnum==3,1:3),'vertices',vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
randColor = rand([1 3]);
patch('faces',tria(tnum==4,1:3),'vertices',vert, ...
    'facecolor',randColor, ...
    'edgecolor',[0,0,0]) ;
title(['MESH-Smoothed.: KIND=DELFRONT, |TRIA|=', ...
    num2str(size(tria,1))]) ;
hold off;

rmpath(modelerPath);
rmpath(mesherPath);