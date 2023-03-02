function opCode = mod2D_opStrToCode(opString)

% Convert to lower case
opString = tolower(opString);

switch opString
  case 'add'
    opCode = 3;
  case 'subtract'
    opCode = 0;
  case 'intersect'
    opCode = 1;
  otherwise
    error('opStrToCode',sprintf('Boolean operation type ''%s'' not supported',opString));
end
