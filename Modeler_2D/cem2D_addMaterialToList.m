function materialList = cem2D_addMaterialToList(material,materialList)
  
  if(~exist('materialList'))
    materialList = [];
  end
  
  nMaterials = numel(materialList);
  
  % If there is no material to add, add the default
  if(~exist('material'))
    material = cem2D_createMaterialDefs;
  end
  
  % If for some reason the user passed an empty material, add the default.
  if(isempty(material))
    material = cem2D_createMaterialDefs;
  end
  
  % Check if material already exists:
  % Extract name
  newMaterialName = material.name;
  % Search for material in list
  existingMaterial = cem2D_getMaterialPropsFromName(newMaterialName,materialList);
  % If exists already, throw an error.
  if(~isempty(existingMaterial))
    errStruct = mod2D_createErrorMsg(...
                  sprintf('Material called %s already exists',newMaterialName),...
                  'cem2D_addMaterialToList');
    throw(errStruct)
  end
  
  % Create a temporary material list
  newMaterialList = cell([nMaterials+1 1]);
  
  % Add all pervious materials
  newMaterialList(1:nMaterials) = materialList;
  
  newMaterialList{end} = material;
  
  materialList = newMaterialList;
  
end