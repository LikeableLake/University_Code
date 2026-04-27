function  dPhidt  = Conv2DFun(~,PHI)

global Nx Ny U V h

dPhidt = zeros(Nx,Ny);
for i=2:Nx-1
    for j=2:Ny-1
        dPhidt(i,j)= - ( U(i,j)*(PHI(i+1,j)-PHI(i-1,j))/(2*h)+...
                         V(i,j)*(PHI(i,j+1)-PHI(i,j-1))/(2*h) );
    end
end

end
