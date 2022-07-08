
% Solves time-dependent lubrication flow over topography 
global m n xvec yvec deltax rho g mu q my mx

% eruption duration (T seconds)
days=5;
T=3600*24*days;

% select pre or post eruption topography
topo_selected=topo_pre;dTdx=dTdx_pre;dTdy=dTdy_pre;
%topo_selected=topo_post;dTdx=dTdx_post;dTdy=dTdy_post;

% parameters
rho=2650;       %density (kg/m3)    }
g=9.81;         %gravity (m/s2)     }-- lambda = rho*g/mu
mu=5e6;         %viscosity (Pas)    }
q=80;           %lava flux (m3/s)  --- input rate = q


% generate rectangular mesh (50m resolution)
deltax= 50; 
xvec=-1700:deltax:1400; n =length(xvec);
yvec=-1300:deltax:1300; m =length(yvec);     
[Xmesh,Ymesh] = meshgrid(xvec,yvec);

% specify the topographic gradient matrices
[xgrid,ygrid]=ndgrid(xdom,ydom);
Tx=griddedInterpolant(xgrid,ygrid,dTdx'); mx=Tx(Xmesh',Ymesh')';
Ty=griddedInterpolant(xgrid,ygrid,dTdy'); my=Ty(Xmesh',Ymesh')';

% initial h=0 everywhere
h00= Xmesh*0;

% timestep using matlab ode solver
options = odeset('NonNegative',[1:m*n]);
[t,hsol]=ode15s(@deriv,[0 T],h00(:),options); 


% plot solution for t=T
hsol_final=reshape(hsol(end,:),[m,n]);
hmax=max(hsol_final,[],'all');
hav=flowaverage(hsol_final,0);
contourh=0.7*hav;

figure; ax=gca; hold on;
% plots lava thickness in red-yellow
contourf(xvec,yvec,hsol_final,linspace(contourh,hmax,100),'Linecolor','none'); % shaded contour plot
colormap autumn; caxis([4 16]); colorbar;
% overlies topographic contours in black
contour(xdom,ydom,topo_selected,[1699:10:2099],'k'); 
% adds observed lava field outline
contour(xdom,ydom,lava_thickness_unsmoothed,[1 1],'linecolor','#32CD32','Linewidth',1.5); %outline of lavafield 

axis([-1600 1600 -1200 1200]);
ax.DataAspectRatioMode='manual';ax.DataAspectRatio = [1 1 1];
title('Lava flow simulation, Marcath');
xlabel('x (m)');
ylabel('y (m)');


% check volume conservation:
err=1-trapz(xvec,trapz(yvec,hsol_final))/(q*T);
err % should be close to zero

%function to discretise spatial derivatives
function dhdt = deriv(t,h)
    global m n xvec yvec deltax rho g mu q mx my
     
    h = reshape(h,[m,n]); % reshape column vector into m x n matrix
    dhdt = zeros(m,n);

    pjplus  = 0.5*( (h(2:m-1,2:n-1).^3).*( mx(2:m-1,2:n-1) + (h(2:m-1,3:n)-h(2:m-1,2:n-1))/deltax ) + (h(2:m-1,3:n).^3).*( mx(2:m-1,3:n) + (h(2:m-1,3:n)-h(2:m-1,2:n-1))/deltax ));
    pjminus = 0.5*( (h(2:m-1,2:n-1).^3).*( mx(2:m-1,2:n-1) + (h(2:m-1,2:n-1)-h(2:m-1,1:n-2))/deltax ) + (h(2:m-1,1:n-2).^3).*( mx(2:m-1,1:n-2) + (h(2:m-1,2:n-1)-h(2:m-1,1:n-2))/deltax ));

    pjplusy  = 0.5*( (h(2:m-1,2:n-1).^3).*( my(2:m-1,2:n-1) + (h(3:m,2:n-1)-h(2:m-1,2:n-1))/deltax ) + (h(3:m,2:n-1).^3).*( my(3:m,2:n-1) + (h(3:m,2:n-1)-h(2:m-1,2:n-1))/deltax ));
    pjminusy = 0.5*( (h(2:m-1,2:n-1).^3).*( my(2:m-1,2:n-1) + (h(2:m-1,2:n-1)-h(1:m-2,2:n-1))/deltax ) + (h(1:m-2,2:n-1).^3).*( my(1:m-2,2:n-1) + (h(2:m-1,2:n-1)-h(1:m-2,2:n-1))/deltax ));

    dhdt(2:m-1,2:n-1) = (rho*g/(3*mu))*((pjplus-pjminus)/deltax + (pjplusy-pjminusy)/deltax) +(q/(3*deltax^2))*((yvec(2:m-1)==0)'*(xvec(2:n-1)==1100)+(yvec(2:m-1)==100)'*(xvec(2:n-1)==1100)+(yvec(2:m-1)==-100)'*(xvec(2:n-1)==1000));

    dhdt=dhdt(:);       %reshapes as column vector with m*n componnents
    
end
