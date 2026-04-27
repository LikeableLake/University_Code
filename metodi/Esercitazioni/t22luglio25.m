clc;clear;close all

% dati

k=1;
AR=2;
lx=1;
ly=AR*lx;
nx=60;
ny=1+AR*(nx-1);
x=linspace(0,lx,nx);
y=linspace(0,ly,ny);
h=x(2)-x(1);
G=1:(nx-2)*(ny-2);
G=reshape(G,ny-2,nx-2);

G=[zeros(ny-2, 1) G zeros(ny-2,1)];
G=[zeros(1,nx); G; zeros(1,nx)];

L=-delsq(G);
hq=h^2;

%vettore termini noti
q=zeros(ny-2, nx-2);
q(:,1)=-k*1/h;
q(1,:)=k*1/hq;
q(end,:)=k*(lx-x(2:nx-1))/hq;

p=q(:);

%%

for i=1:ny
    for j=2:nx-1
        node=G(i,j);
        if (G(i, j-1)==0 && G(i, j)~=0)
            L(node,node)=-3;
        end
        if (G(i, j+1)==0 && G(i, j)~=0)
            L(node, node)=-3;
        end
    end
end


%soluzione


PHI=-(L\p)*hq/k;

PHI=reshape(PHI, ny-2,nx-2);

PHI=[PHI(:,1), PHI, PHI(:,end)];
PHI=[(ones(1, nx)); PHI; (lx-x)];

figure(1);
surf(x,y,PHI, 'LineStyle', 'none', 'FaceColor', 'interp'); 
xlabel('x');
ylabel('y');
zlabel('\Phi ', 'Rotation', 0);
colorbar;


%%

beta=0.8;
k=1;
dt=hq*beta/k;
T=(lx^2)*beta/k;
nt=round(T/dt);
dt=T/nt;

PHI0=zeros(ny,nx);

for i= 1:ny
    for j= 1:nx
        PHI0(i , j)= y(i)*(lx-x(j))/(2*lx);
    end
end

%% (I-beta/2 * LAP)PHIn+1=(I+beta/2 * LAP) PHIn + P

I=eye(size(L), "like", L);

A=I+((beta/2)*L);
B=I-((beta/2)*L);

%% condizioni al contorno
q=zeros(ny-2, nx-2);
q(:,1)=-k*1*h*beta;
q(1,:)=beta*1;
q(end,:)=beta*(lx-x(2:nx-1));

p=q(:);

%%

%BPHIn+1=APHIn+P

PHIf=PHI0;
PHIf=PHIf(2:ny-1, 2:nx-1);
PHIf=PHIf(:);

for it=1:400
    

    PHIf=B\(A*PHIf+p);


    Phif=reshape(PHIf, ny-2, nx-2);
    Phif=[Phif(:,1), Phif, Phif(:,end)];
    Phif=[(ones(1, nx)); Phif; (lx-x)];

    figure(2);
    surf(x,y,Phif, 'LineStyle', 'none', 'FaceColor', 'interp'); 
    xlabel('x');
    ylabel('y');
    zlabel('\Phi ', 'Rotation', 0);
    colorbar; 
    drawnow;
    hold off;
end

   
