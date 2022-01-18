function theta=ori_map(Nx,Kc)
% Formula from Kaschube et al. (2010) Supp Materials Eq. 20. 
% Kaschube, M., Schnabel, M., L¨owel, S., Coppola, D. M., White, L. E., and Wolf, F. (2010).
% Universality in the evolution of orientation columns in the visual cortex. Science, 330:1113?
% 1116.

% Nx=50;
Ny=Nx;
dx=1/Nx; 
x=0:dx:1;
y=x;

n=30; 
% Kc=5*2*pi; % spatial frequency 
k=Kc*[cos((0:n-1)'*pi/n), sin((0:n-1)'*pi/n)];
l=sign(rand(n,1)-.5); l(l==0)=1; 
phi=rand(n,1)*2*pi; 
z=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        z(i,j)=sum(exp(sqrt(-1)*(l.*(k*[x(i); y(j)])+phi))); 
    end
end
theta=angle(z)/pi/2+.5; % preferred orientation, theta in [0 1]
% figure
% imagesc(theta); 
% colormap(hsv)