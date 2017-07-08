function materialStruct = cem2D_getMaterialPropsFromName(materialName,materialList)

% In order for the 'for' loop to work, this must be a column cell vector
if(size(materialList,2) ~= 1)
  materialList = materialList(:);
end

materialStruct = [];

for cMaterialCell = materialList;
  cMaterial = cMaterialCell{1};
  % Check if this is 
  if(strcmp(cMaterial.name,materialName))
    materialStruct = cMaterial;
    break;
  end
end

% If there is no material with this name, then break it up

end