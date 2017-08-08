function materialStruct = cem2D_getMaterialPropsFromName(materialName,materialList)

% Check if there are no materials
if(isempty(materialList))
  materialStruct = [];
  return;
end

% In order for the 'for' loop to work, this must be a column cell vector
if(size(materialList,2) ~= 1)
  materialList = materialList(:);
end

materialStruct = [];


for materialIdx = 1:numel(materialList)
  cMaterial = materialList{materialIdx};
  % Check if this is 
  if(strcmp(cMaterial.name,materialName))
    materialStruct = cMaterial;
    break;
  end
end

% If there is no material with this name, then break it up

end