% test_matrix_construction_material

clc
clear
close all

%% Load relevant data

modelerPath = './Modeler_2D';
mesherPath = '/home/gadi/Repositories/mesh2d';
meshWrapperPath = './Mesh_2D';
solverPath = './Solver_2D';

addpath(modelerPath);
addpath(mesherPath);
addpath(meshWrapperPath);
addpath(solverPath);

load solver_predata.mat

%% Construct matrix for volume\material B.Cs.

[K,b] = cem2D_createKmatBvect_materials(...
          smoothMesh,...            % Mesh data structure
          materialList,...          % Cell array of all materials
          materialAssignment,...    % Cell array assigning materials to faces
          simProps, ...             % Simulation properties
          1.5);                     % Current simulation frequency (in project units)



%%

%rmpath(modelerPath);
%rmpath(mesherPath);
%rmpath(meshWrapperPath);
%rmpath(solverPath);