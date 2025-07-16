function Uxx = Laplace_solver(x,y,dx,dy,f,duL,uR,uD,uT)

format long

Nx = length(x); Ny = length(y);
% Create matrices of Mu=F
M = zeros(Nx*Ny); F = zeros(Nx*Ny,1);       

%Solving interior
for i = 2:Nx-1      
    for j = 2:Ny-1
        k=(j-1)*Nx+i;
        M(k,[k+Nx k k-Nx])=M(k,[k+Nx k k-Nx]) + [1 -2 1]/dy^2;
        M(k,[k+1 k k-1])=M(k,[k+1 k k-1]) + [1 -2 1]/dx^2;
        F(k)=f(i,j);
    end
end
%Solving Left Boundary
i=1; 
for j=1:Ny     
    k=(j-1)*Nx+i;
    M(k,[k+1 k])=M(k,[k+1 k]) +[1 -1]/(dx);
    F(k)= duL(j);
end
%Solving Right Boundary
i=Nx; 
for j=1:Ny     
    k=(j-1)*Nx+i;
    M(k,[k])=M(k,[k]) + [1];
    F(k)=uR(j);
end
%Solving Bottom Boundary
j=1; 
for i=1:Nx     
    k=(j-1)*Nx+i;
    M(k,[k])=M(k,[k]) + [1];
    F(k)=uD(i);
end
%Solving Top Boundary
j=Ny;
for i=1:Nx     
    k=(j-1)*Nx+i;
    M(k,[k])=M(k,[k]) + [1];
    F(k)=uT(i);
end
% Solvng Mu = F 
Uxx = M\F;

end

