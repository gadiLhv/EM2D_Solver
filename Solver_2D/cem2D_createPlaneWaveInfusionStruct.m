function pwStruct = cem2D_createPlaneWaveInfusionStruct(varargin)
% Polarization is set via simulation properties

    pwStruct = struct('phiPropagation',0,...  % Phi of propagation direction, in degrees
                      'amp',1,...             % Amplitude in V/m (for E-field. Need to be corrected for 0 field
                      'polarizationVect',[1 0],...  % Computed around propagatin axis, in transverse to the computation plane
                      'measurementPlaneDist',0,...  % Phase is zero at this plane. Defaulted in [0,0].
                      'measurementStartingPoint',[0 0],... % Distance is measured from this coordinate
                      );
        

    pwStruct = misc_validatePropStruct(pwStruct,varargin);

  
  
end
