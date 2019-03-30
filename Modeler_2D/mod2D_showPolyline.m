function polHdl = mod2D_showPolyline(axHdl,lineStruct,edgeColor)
  hold(axHdl)
  
  % Set right the rotation directions (CCW, fill in, CW, fill out)
  plot (lineStruct.x,lineStruct.y,...
         '-','linewidth',3,...
         'color', edgeColor);
  
  hold off;
end