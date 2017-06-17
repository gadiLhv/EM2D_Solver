function normFactor = mod2D_lengthUnitFactor(unitTxt)
  unitTxt = validatestring( unitTxt,...
                            {'nm','um','mm','cm','m'});
                            
  switch(unitTxt)
    case 'nm'
      normFactor = 1e9;
    case 'um'
      normFactor = 1e6;
    case 'mm'
      normFactor = 1e3;
    case 'cm'
      normFactor = 1e2;
    case 'm'
      normFactor = 1;
  end
end
