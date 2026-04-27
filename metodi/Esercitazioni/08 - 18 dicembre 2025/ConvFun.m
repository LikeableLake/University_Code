function dPhidt = ConvFun(~,Phi)

global A

Phi    = Phi(:);
dPhidt = A*Phi;

end
