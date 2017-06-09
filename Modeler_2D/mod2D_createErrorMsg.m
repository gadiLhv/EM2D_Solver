function errMsgStruct = mod2D_createErrorMsg(errMsg,invocation)
  errMsgStruct = struct('msg',errMsg,'caller',invocation);
end