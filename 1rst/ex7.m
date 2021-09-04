% ------------------------------------------------
% Define computational parameters
% ------------------------------------------------
 a=0.1;
 d=1; 
 h=1;
 nfig=1; 
 e0 = 8.85e-12;
 Dx=0.01;
 Dy=0.01;
[XP,Y] = meshgrid(linspace(0,2,100));
[XS,Z] = meshgrid(linspace(0,2,100),linspace(-2,2,100));
% ----------------------------------------------------
% array calculation for potential 
% (faster) numerical integration
% ----------------------------------------------------
h=1;
Phi_1=a*integral(@(phi)1./R_vector(XP,Y,d,h,a,phi),0,2*pi,'ArrayValued',true);
Phi_2=-a*integral(@(phi)1./R_vector(XP,Y,d,-h,a,phi),0,2*pi,'ArrayValued',true);
Phi_3=a*integral(@(phi)1./R_vector(XP,Y,-d,-h,a,phi),0,2*pi,'ArrayValued',true);
Phi_4=-a*integral(@(phi)1./R_vector(XP,Y,-d,h,a,phi),0,2*pi,'ArrayValued',true);
Phi = Phi_1+Phi_2+Phi_3+Phi_4;
Phi_max = max(max(Phi));
% -----------------------------
% Plot Potential
% -----------------------------
figure(nfig)
surface(XP,Y,Phi,'LineStyle','None');
hold on 
set(gca,'Fontsize',12,'Fontweight','bold')
xlabel('x','Fontsize',12,'FontWeight','bold')
ylabel('y','Fontsize',12,'FontWeight','bold')
title(['Normalized Potential \Phi (x,y)/(\lambda/4\pi\epsilon_0)',' = 1C/m, d = �
1m, h = 1m, a = 0.1m'], 
 'Fontsize',10,'FontWeight','bold','Color','b')
 
hold on 
t1=linspace(0.9,1.1,100);
y9=1; 
w=plot(t1,y9,'k','Linewidth',2,'Color','r');
hold on
nfig=nfig+1; 
figure(nfig) 
cont = [0.01 0.1 0.2 0.4 0.8 1.0 1.5 2 2.5 3 3.5 4 5 6 7.5];
[CS,H] = contour(XP,Y,Phi,cont,'Linewidth',1); hold on 
clabel(CS,H,cont);
hold on 
set(gca,'Fontsize',12,'Fontweight','bold')
xlabel('x','Fontsize',12,'FontWeight','bold')
ylabel('y','Fontsize',12,'FontWeight','bold')
title(['Normalized Potential \Phi (x,y)/(\lambda/4\pi\epsilon_0)',' = 1C/m, d = �
1m, h = 1m, a = 0.1m'], 
 'Fontsize',10,'FontWeight','bold','Color','b')
axis equal
% ----------------------------------------------------
% Calculation of electric field
% using numerical differention (gradient)
% ----------------------------------------------------
[Exarray1,Eyarray1]=gradient(Phi,Dx,Dy);
Exarray = - Exarray1;
Eyarray = - Eyarray1;
% --------------------------------------------------
% Plot Electric Field (quiver & streamslice)
% --------------------------------------------------
nfig = nfig + 1;
figure(nfig)
hold off
LL = ((Exarray).^2 + (Eyarray).^2).^(1/2);
quiver(XP,Y,Exarray./LL,Eyarray./LL,0.5) 
hold on 
t1=linspace(0.9,1.1,100);
y9=1; 
plot(t1,y9,'k','Linewidth',2,'Color','r');
hold on
hs = streamslice(XP,Y,Exarray,Eyarray,2);
set(hs,'Color','m','Linewidth',1.0)
contour(XP,Y,Phi,cont,'Linewidth',1,'Color','b');
hold on 
t1=linspace(0.9,1.1,100);
y9=1; 
plot(t1,y9,'k','Linewidth',2,'Color','r');
hold on
set(gca,'Fontsize',12,'Fontweight','bold')
xlabel('x','Fontsize',12,'FontWeight','bold')
ylabel('y','Fontsize',12,'FontWeight','bold')
title(['Normalized Electric Field E_x,E_y /(\lambda/4\pi\epsilon_0)',' = 1C/m, d �
= 1m, h = 1m, a = 0.1m'], ...
 'Fontsize',10,'FontWeight','bold','Color','b')
axis equal
grid on
% --------------------------------------------------
% array calculation for surface 
% --------------------------------------------------
h=1;
factor = h*e0;
sur_1=-factor.*integral(@(phi)1./
Rn_vector(XS,Z,d,h,a,phi),0,2*pi,'ArrayValued',true);
sur_2=-factor.*integral(@(phi)1./Rn_vector(XS,Z,d,-
h,a,phi),0,2*pi,'ArrayValued',true);
sur_3=factor.*integral(@(phi)1./Rn_vector(XS,Z,-d,-
h,a,phi),0,2*pi,'ArrayValued',true);
sur_4=factor.*integral(@(phi)1./Rn_vector(XS,Z,-
d,h,a,phi),0,2*pi,'ArrayValued',true);
sur=sur_1+sur_2+sur_3+sur_4;
nfig=nfig+1;
figure(nfig) 
contour(XS,Z,sur,'Showtext','on','Linewidth',1);
set(gca,'Fontsize',12,'Fontweight','bold')
xlabel('x','Fontsize',12,'FontWeight','bold')
ylabel('z','Fontsize',12,'FontWeight','bold')
title(['Normalized \ (x,y)/(\ *a/4\pi\epsilon_0)',' = 1C/m, d = 1m, h = 1m, a = � � �
0.1m'], 
 'Fontsize',10,'FontWeight','bold','Color','b')
axis equal 
hold on 
t = linspace(0,2*pi,100); 
circsy2 = a.*sin(t)+1; 
circsz2 = a.*cos(t); 
p1 = plot(circsy2,circsz2,'k','Linewidth',2,'Color','r'); 
% =====================================================================
function R = R_vector(x,y,d,h,a,phi)
 R = (((x-d).^2+(y-h).^2+a^2-2*a.*(x-d).*cos(phi)).^(1/2);
end
% =====================================================================
function Rn = Rn_vector(x,z,d,h,a,phi)
 Rn = ((x-d).^2+h^2+a^2+z.^2 - 2*a.*((x-d).*cos(phi)+z.*sin(phi))).^(3/2);
end
% ====================================================================