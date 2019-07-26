function [materialName,materialProps] = misc_chooseWeakestMaterial(material1name,material1defs,material2name,material2defs)
    materialRating = {'PEC','PMC','metal','normal'};
    get1Rating = @(matType) strcmp(material1defs.type,matType);
    get2Rating = @(matType) strcmp(material2defs.type,matType);
    
    rating1 = find(cellfun(get1Rating,materialRating));
    rating2 = find(cellfun(get2Rating,materialRating));
    
    if rating1 < rating2
        materialName = material2name;
        materialProps = material2defs;
    else
        materialName = material1name;
        materialProps = material1defs;
    end
end
