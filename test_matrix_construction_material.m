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
% Change polarization, for the heck of it
simProps.polarizationType = 'TE';
%% Construct matrix for volume\material B.Cs.

[K,b] = cem2D_createKmatBvect_materials(...
          smoothMesh,...            % Mesh data structure
          materialList,...          % Cell array of all materials
          materialAssignment,...    % Cell array assigning materials to faces
          simProps, ...             % Simulation properties
          1.5);                     % Current simulation frequency (in project units)

[K,b] = cem2D_createKmatBvect_2ndOrderRadCond(...
          smoothMesh,...            % Mesh data structure
          materialList,...          % Cell array of all materials  
          materialAssignment,...    % Cell array assigning materials to faces
          simProps, ...             % Simulation properties
          1.5);                     % Current simulation frequency (in project units)


pwStruct = struct('phiPropagation',0,...  % Phi of propagation direction
        'amp',1,...             % Amplitude in V/m.
        'polarizationVect',[1 0],...  % Computed around propagatin axis, in transverse to the computation plane
        'measurementPl aneDist',0,...  % Phase is zero at this plane. Defaulted in [0,0].
        'measurementStartingPoint',[0 0],... % Distance is measured from this coordinate
        );
          
pwProps = cem2D_createPlaneWaveInfusionStruct(...
            'phiPropagation',45)
cem2D_createPlaneWaveInfusionStruct
          
%%

rmpath(modelerPath);
rmpath(mesherPath);
rmpath(meshWrapperPath);
rmpath(solverPath);