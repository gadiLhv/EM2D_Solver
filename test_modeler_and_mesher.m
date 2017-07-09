% test_modeler_and_mesher
clc
clear
close all

modelerPath = './Modeler_2D';
mesherPath = '/home/gadi/Documents/Octave_Projects/Mesh2D_Latest';

%% Mesher parameters
% Set printout/no printout
outputMesh = true(1);

% Global mesh properties

f0 = 2e9;
c0 = 3e8;
% 1. Maximum allowed mesh size with respect to wavelength.
% Take note that this takes into consideration the material
% properties.
relWLmeshMax = 0.33;

% 2. Bounding box addition (Currently this is mesh
% properties, in the future this needs to belong to B.C.)
% This is given in relaitve terms of wavelength.
boundingBoxAddSpace = 0.25;

% Create default material


addpath(modelerPath);
addpath(mesherPath);


% Create the default material
defaultMaterial = cem2D_createMaterialDefs;
% Create the material list with only the default material
materialList = cem2D_addMaterialToList(defaultMaterial);

% Add three materials
material3p1 = cem2D_createMaterialDefs('name','plastic310','er',3.1);
material5p5 = cem2D_createMaterialDefs('name','plastic550','er',5.5);
material7p2 = cem2D_createMaterialDefs('name','plastic720','er',7.2);

materialList = cem2D_addMaterialToList(material3p1,materialList);
materialList = cem2D_addMaterialToList(material5p5,materialList);
materialList = cem2D_addMaterialToList(material7p2,materialList);

%% Create the geometry


% Geometry 1
pol1 = mod2D_createRectangleStruct([-2 -2]*10,[2 -1]*10);
pol2 = mod2D_createRectangleStruct([-2 -1]*10,[2 0]*10);
pol3 = mod2D_createRectangleStruct([-2 0]*10,[2 1]*10);
rect1 = mod2D_createRectangleStruct([-1 -1.5]*10,[1 0.5]*10);
pol1 =  mod2D_booleanOperation(pol1,rect1,'subtract');
pol2 =  mod2D_booleanOperation(pol2,rect1,'subtract');
pol3 =  mod2D_booleanOperation(pol3,rect1,'subtract');

% Geometery 2
%rect1 = mod2D_createRectangleStruct([-3 -3],[3 3]);
%rect2 = mod2D_createRectangleStruct([-2.5 -2.5],[2.5 2.5]);
%pol1 = mod2D_booleanOperation(rect1,rect2,'subtract');
%rect1 = mod2D_createRectangleStruct([-2 -2],[2 2]);
%rect2 = mod2D_createRectangleStruct([-1 -1],[1 1]);
%pol2 = mod2D+_booleanOperation(rect1,rect2,'subtract');
%pol1 = mod2D_booleanOperation(pol1,pol2,'add');
%rect1 = mod2D_createRectangleStruct([2 -0.5],[2.5 1.5]);
%pol1 = mod2D_booleanOperation(pol1,rect1,'add');
%pol2 = mod2D_createRectangleStruct([-2.5 1],[-2 2]); % Bad at handling duplicate nodes. Remove them while creating connectivty list
%pol3 = mod2D_createPolygonStruct([-0.7 0.4 1],[-0.6 -0.8 1]);

% Create polygon list
polList = {pol1,pol2,pol3};

% Get default material properties
defaultMaterial = cem2D_getMaterialPropsFromName('default',materialList);
default_er = defaultMaterial.er;
default_mr = defaultMaterial.mr;

% Create bounding box (background)
bBox = mod2D_createBoundingBox(...
        polList,...
        boundingBoxAddSpace*(c0/f0)/sqrt(default_er*default_mr));

% Assign materials to polygons


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

% Diagonal of the bounding box is the maximum possible distance
maxEdgeL = sqrt((bgBoundX(2)-bgBoundX(1)).^2 + (bgBoundY(2)-bgBoundY(1)).^2);


% Make sure that no new "free" nodes are added. In case that new nodes are added
% in "mid-air", they need to be determined inside the polygon. 
lsfOpts.kind = 'delaunay';
lsfOpts.dhdx = 1;
lsfOpts.rho2 = sqrt(maxEdgeL*10);
% Estimate Local Feature Size (LFS) for each node\part
[vlfs,tlfs,hlfs] = lfshfn2( node,...
                            edge,...
                            face,...
                            lsfOpts);

           
           
% Show initial feature size estimation
patch('faces',tlfs,'vertices',vlfs,'facecolor',[1 1 1],'edgecolor',[0 0 0])
hold on;
plot(node(:,1),node(:,2),'.r','markersize',15);
hold off;

%hs = [sqrt(sum((vlfs(tlfs(:,2),:)-vlfs(tlfs(:,1),:)).^2,2)) ...
%      sqrt(sum((vlfs(tlfs(:,3),:)-vlfs(tlfs(:,2),:)).^2,2)) ... 
%      sqrt(sum((vlfs(tlfs(:,1),:)-vlfs(tlfs(:,3),:)).^2,2))];


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
%[ vert,...              % Vertices of mesh cells
%  etri,...              % Constrained edge list (???)
%  tria,...              % Triangle threesomes (attached to VERT)
%  tnum] = ...           % Part assignments
%  refine2(node, ...     % Full node list
%          edge, ...     % Edge connectivity
%          face);

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