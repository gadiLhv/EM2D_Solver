function materialList = cem2D_addMaterialToList(material,materialList)
  
  if(~exist('materialList'))
    materialList = [];
  end
  
  nMaterials = numel(materialList);
  
  
  
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
  
  errMsgStruct = mod2D_createErrorMsg(errMsg,invocation)
  % Create a temporary material list
  newMaterialList = cell([nMaterials+1 1]);
  
  % Add all pervious materials
  newMaterialList(1:nMaterials) = materialList;
  
  newMaterialList{end} = material;
  
  materialList = newMaterialList;
  
end