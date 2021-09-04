function Exercise_6(a,b,d,D,L,nfig)
% ------------------------------------
% Exercise_6(10,5,3,2.5,3,1)
% ------------------------------------
% ------------------------------------------------
% Define computational parameters
% ------------------------------------------------
%a=10;
%b=5;
%d=3;
%L=3; 
%nfig=1;
epsilon0 = 8.854e-12;
ymin=-15;
ymax=15; 
zmin=-15;
zmax=15;
Dy=(ymax-ymin)/200;
Dz=(zmax-zmin)/200;
y1=ymin:(ymax-ymin)/200:ymax; 
z1=zmin:(zmax-zmin)/200:zmax; 
[Y,Z] = meshgrid(y1,z1);
Farray = zeros(length(y1),length(z1));
Theta= L/d; 
% ----------------------------------------------------
% Array calculation for potential 
% numerical integration
% ----------------------------------------------------
%Calculation of the Potential 
for iy = 1:length(y1)
 for iz = 1:length(z1)
 yy1=Y(iy,iz);
 zz1=Z(iy,iz);
 [func] = potential(yy1,zz1,a,d,D,b,Theta);
 Farray(iy,iz) = func; 
 end 
end 
Phi_max = max(max(Farray));
% -----------------------------
% Plot Potential
% -----------------------------
figure(nfig) 
hold off 
surface(Y,Z,Farray), shading interp
hold on
t = linspace(0,2*pi,100); 
circsy1 = a.*sin(t); 
circsz1 = a.*cos(t); 
p1 = plot(circsy1,circsz1,'k','Linewidth',2); 
circsy2 = b.*sin(t); 
circsz2 = b.*cos(t) + D; 
p2 = plot(circsy2,circsz2,'k','Linewidth',2); 
tarc = -Theta:2*Theta/100:Theta;
yarc = d*sin(tarc);
zarc = d*cos(tarc) +D;
p3 = plot(yarc,zarc,'r','Linewidth',2); 
set(p1,'ZData',1.5*Phi_max+1+zeros(size(circsy1)));
set(p2,'ZData',1.5*Phi_max+1+zeros(size(circsy2)));
set(p3,'ZData',1.5*Phi_max+1+zeros(size(yarc)));
caxis([0 6]);
colorbar
axis equal
cont =[0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.10 1.25 2 3 5];
nfig=nfig+1;
figure(nfig);
hold off
[CS,H] = contour(Y,Z,Farray,cont,'Linewidth',1);
clabel(CS,H,cont);
axis equal
% ----------------------------------------------------
% Calculation of electric field
% using numerical differention (gradient)
% ----------------------------------------------------
[Exarray1, Eyarray1] = gradient(Farray,Dy,Dz);
Exarray = - Exarray1;
Eyarray = - Eyarray1;
% --------------------------------------------------
% Plot Electric Field (quiver & streamslice)
% --------------------------------------------------
nfig = nfig + 1;
figure(nfig)
hold off
LL = ((Exarray).^2 + (Eyarray).^2).^(1/2);
quiver(Y,Z,Exarray./LL,Eyarray./LL,0.5); 
hold on
hs = streamslice(Y,Z,Exarray,Eyarray,2);
set(hs,'Color','m','Linewidth',1.0);
hold on
contour(Y,Z,Farray,cont,'Linewidth',1,'Color','b');
axis equal
% ----------------------------------------------------------------
% Array Calculation of electric field
% for the surface density on the verge of the two surfaces 
% ----------------------------------------------------------------
rho=5; 
phi=0:2*pi/200:2*pi; 
z8=rho*cos(phi)+D;
y8=rho*sin(phi);
 Eyarray = integral(@(theta) electrical_field_y(theta,y8,z8,d,D,b),-
Theta,Theta,'ArrayValued',true,...
 'RelTol',1e-6,'AbsTol',1e-12);
 Ezarray = integral(@(theta) electrical_field_z(theta,y8,z8,d,D,b),-
Theta,Theta,'ArrayValued',true,...
 'RelTol',1e-6,'AbsTol',1e-12);
% --------------------------------------------------
% Array calculation for Surface Density
% --------------------------------------------------
nfig=nfig+1;
figure(nfig)
[s_c] = surface_charge(phi,b,D,Eyarray,Ezarray);
plot(phi*180/pi,s_c,'Linewidth',2)
set(gca,'Fontsize',12,'Fontweight','bold')
xlabel('\theta (deg)','Fontsize',12,'FontWeight','bold')
ylabel('Normalized Surface Charge, \sigma_b/\lambda 
(1/m)','Fontsize',12,'FontWeight','bold')
grid on
end
% ==========================================================================
% --------------------------------------------------
% Functions 
% --------------------------------------------------
function [func] = potential(y5,z5,a,d,D,b,Theta)
Rc = ((z5-D)^2+y5^2)^(1/2);
r = (z5^2 + y5^2)^(1/2);
if (Rc<=b)
 func = integral(@(theta) potential_1(theta,y5,z5,d,D,b),-
Theta,Theta,'RelTol',1e-6,'AbsTol',1e-12);
 func = func+1;
end
if (Rc>b && r<=a)
 func=1;
end
if (r>a)
 func=a/r;
end
endfunction
% ===========================================================================
function [dphi]=potential_1(theta,y,z,d,D,b)
ds=b^2/d;
R = ((y-d*sin(theta)).^2+(z-d*cos(theta)-D).^2).^(1/2);
R1= ((y-ds*sin(theta)).^2 + (z-ds*cos(theta)-D).^2).^(1/2);
dphi=d./R - b./R1; % + 1; 
endfunction 
% ===========================================================================
function [delec]=electrical_field_y(theta,y,z,d,D,b)
ds=b^2/d;
R = ((y- d*sin(theta)).^2 + (z- d*cos(theta)-D).^2).^(1/2);
R1= ((y-ds*sin(theta)).^2 + (z-ds*cos(theta)-D).^2).^(1/2);
rn= (y-d*sin(theta));
Rn= (y-ds*sin(theta));
delec=d*rn./R.^3 - b*Rn./R1.^3; 
endfunction
% ==========================================================================
function [delec]=electrical_field_z(theta,y,z,d,D,b)
ds=b^2/d;
R = ((y- d*sin(theta)).^2 + (z- d*cos(theta)-D).^2).^(1/2);
R1= ((y-ds*sin(theta)).^2 + (z-ds*cos(theta)-D).^2).^(1/2);
rn= (z-d*cos(theta)-D);
Rn= (z-ds*cos(theta)-D);
delec=d*rn./R.^3 - b*Rn./R1.^3; 
endfunction
% ========================================================================== 
function [s_c]=surface_charge(theta,b,D,Eyarray,Ezarray)
y=b*sin(theta);
z=b*cos(theta)+D;
rb=(y.^2+(z-D).^2).^(1/2);
s_c=-((y./rb).*Eyarray + ((z-D)./rb).*Ezarray);
s_c=s_c/(4*pi);
endfunction
% ========================================================================== 
