% solves lubrication flow over topography


syms m x y u(x,y) D L tanb

% parameters describing topography - change these to explore the effect of
% changes in topography
M=0.5;      % ratio of the aspect ratio of the feature to the background slope
beta=10;    % background slope in degrees
L=100;      % feature radius - note integration domain will need to change size if this is increased

tanb=tan(beta*pi/180);
D=M*L*tanb;

% symbolic description of topography (see section 3.2 in manuscript)
m(x,y)=-x.*tanb+D*(x.^2+y.^2).*exp((-x.^2-y.^2)./L^2)./L^2;  

ms=string(diff(u.^3.*diff(m,x),x)+diff(u.^3.*diff(m,y),y));
ms=strrep(strrep(strrep(strrep(strrep(ms,'^','.^'),'*','.*'),'diff(u(x, y), y)','uy'),'diff(u(x, y), x)','ux'),'u(x, y)','u');

% define geometry for the solver
R1 = [3,4,-300,-300,1050,1050,0,450,450,0]';
geom = R1;
ns = (char('R1'))';
sf = 'R1';
gd = decsg(geom,sf,ns);  

% upstream depth - change to explore the effect of changes in flux, lava
% density and viscosity 
hin=5;

% small additional flux - only required for M>M_2
eps = 0e-5;

% meshing
[p,e,t] = initmesh(gd);        
[p,e,t] = refinemesh(gd,p,e,t);
[p,e,t] = refinemesh(gd,p,e,t);
[p] = jigglemesh(p,e,t);

% max number of triangles for initial
Mt=10000;

% solves −∇⋅(c∇u) + au = f with boundary conditions b in @pdebound
% flow thickness is iteratively decreased n times by a factor of k to approach Hinf
k=1.5; n=3;
global Hinf; 
% initial guess
u0=(hin*k^n).*ones(1,length(p(1,:)));
for i = 1:n
    
    disp(i)
    Hinf=hin*k^(n-i)

    b=@pdebound;
    f=strcat('(',ms,')+',num2str(eps),'.*(heaviside(x+',num2str(L),').*heaviside(',num2str(L),'-x).*heaviside(',num2str(L),'-y))/u');
    a=0;
    c='(u.^3)';

    [u,p,e,t]=adaptmesh(gd,b,c,a,f,'Nonlin','on','Init',u0,...
    'Mesh',p,e,t,'MesherVersion','R2013a','Toln',1e-4,'Ngen',30,'Maxt',...
    Mt,'Jac','full','tripick','pdeadworst','Par',0.5);

    u0=u./k;            

end

% increase resolution of mesh
while length(u)<1e5
    [p,e,t,u]=refinemesh(gd,p,e,t,u);
end
u0=u;

% single iteration of solver
global Hinf; Hinf=hin;

b=@pdebound;
f=strcat('(',ms,')+',num2str(eps),'.*(heaviside(x+',num2str(L),').*heaviside(',num2str(L),'-x).*heaviside(',num2str(L),'-y))/u');
a=0;
c='(u.^3)';

[u,p,e,t]=adaptmesh(gd,b,c,a,f,'Nonlin','on','Init',u0,...
'Mesh',p,e,t,'MesherVersion','R2013a','Toln',1e-4,'Ngen',30,'Maxt',...
Mt,'Jac','full','tripick','pdeadworst','Par',0.5);

% plot as contour plot
figure;
pdeplot(p,e,t,'xydata',u,'contour','on','colorbar','on');
colormap autumn
axis([-300 1000 0 400]);
xlabel('x (m)')
ylabel('y (m)')

% add in outline of topographic feature
th=linspace(0,pi,1000);
xc=L.*cos(th);
yc=L.*sin(th);
hold on; plot(xc,yc,'w','LineWidth',1)
