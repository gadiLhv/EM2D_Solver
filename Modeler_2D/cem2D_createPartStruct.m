function partStruct = cem2D_createPartStruct(polStuct,partData)

% Basic checks
if(isempty(polStruct))
  error('No polygon defined for part');
end

% If there is no supplied part data, use default data
if(isempty(partData))
  partData = CEM2D_createPartDataStruct;
end

partStruct = struct('polygon',polStruct,...
                    'data',partData);

end