%function askhsh_8(h,a,nfig)
%function askhsh_8(1,0.1,0)

%Constants 

d=2;
h = 1; 
a = 0.1; 
nfig = 0;
I=1;
mu_0 = 4*pi*10^(-7);
m0 = I*pi*a^2;

%Create Grid dimensions 
xmin = 0;
xmax = 10;
ymin = 0;
ymax = 10;
zmin = 0;
zmax = 10;

Npoints = 250;
xx = xmin:(xmax-xmin)/Npoints:xmax;
yy = ymin:(ymax-ymin)/Npoints:ymax;
zz = zmin:(zmax-zmin)/Npoints:zmax;

[X, Y] = meshgrid(xx, yy);
[X1, Z] = meshgrid(xx, zz);
[Y1,Z1] = meshgrid(yy,zz);


%Calcualtion of Vector Potential and depiction through quiver and
%streamslice 
Ax = Cal_Of_Vec_Pot_x(X,Y,2,d,h,mu_0,m0);
Az =  Calc_of_Vec_Pot_z(X,Y,2,d,h,mu_0,m0);

nfig = nfig + 1;
figure(nfig)
LL = (Ax.^2 + Az.^2).^(1/2); %kanonikopoihsh 
quiver(X, Y, Ax./LL, Az./LL, 0.5);
ylabel('y axis')
xlabel('x axis')
title('Vector_Potential A(x,y) with quiver in xy plane')
axis equal
grid on
hold off


nfig = nfig + 1;
figure(nfig)
hs=streamslice(X, Y, Ax, Az, 2); 
set(hs,'Color','m','Linewidth',0.8)
ylabel('y axis')
xlabel('x axis')
title('Vector Potential A(x,y) with streamslice in xy plane')   
axis equal
grid on
hold off


Ax1 =  Cal_Of_Vec_Pot_x(X1,1,Z,d,h,mu_0,m0);
Az1 =  Calc_of_Vec_Pot_z(X1,1,Z,d,h,mu_0,m0);


nfig = nfig + 1;
figure(nfig)
LL = (Ax1.^2 + Az1.^2).^(1/2);
quiver(X1, Z, Ax1./LL, Az1./LL, 0.5);
ylabel('z axis')
xlabel('x axis')
title('Vector_Potential A(x,z) with quiver in xz plane')
axis equal
grid on
hold off

 
nfig = nfig + 1;
figure(nfig)
hs=streamslice(X1, Z, Ax1, Az1, 2); 
set(hs,'Color','m','Linewidth',0.8)
ylabel('z axis')
xlabel('x axis')
title('Vector_Potential A(x,z) with streamslice in xz plane')
axis equal
grid on
hold off

Hx =  Calc_Magnetic_Field_x(X,Y,d,h,m0);
Hy =  Calc_Magnetic_Field_y(X,Y,d,h,m0);


nfig = nfig + 1;
figure(nfig)
LL = (Hx.^2 + Hy.^2).^(1/2);%kanonoikopoihsh 
quiver(X, Y, Hx./LL, Hy./LL, 0.5);
ylabel('y axis')
xlabel('x axis')
title('Magnetic Field H(x,y)with quiver in xy plane')
axis equal
grid on
hold off

nfig = nfig + 1;
figure(nfig)
hs=streamslice(X, Y, Hx, Hy, 2); 
set(hs,'Color','m','Linewidth',0.8)
ylabel('y axis')
xlabel('x axis')
title('Magnetic Field H(x,y)with quiver in xy plane')
axis equal
grid on
hold off

Ky =  Calc_Surface_Current_y(0,Y1,Z1,d,h,m0);
Kz =  Calc_Surface_Current_z(0,Y1,Z1,d,h,m0);


nfig = nfig + 1;
figure(nfig)
LL = (Ky.^2 + Kz.^2).^(1/2);
hs=quiver(Y1, Z1, Ky./LL, Kz./LL, 0.5);  
set(hs,'AutoScale', 'on','AutoScaleFactor', 2);
ylabel('y axis')
xlabel('x axis')
title('Surface Current K(y,z) with quiver')
axis equal square
grid on
hold off


nfig = nfig + 1;
figure(nfig)
hs=streamslice(Y1, Z1, Ky, Kz, 2); 
set(hs,'Color','m','Linewidth',0.8)
ylabel('y axis')
xlabel('x axis')
title('Surface Current K(y,z) with streamslice')
axis equal
grid on
hold off

Kx2 =  Calc_Surface_Current_new_x(X1,0,Z,d,h,m0);
Kz2 =  Calc_Surface_Current_new_z(X1,0,Z,d,h,m0);

nfig = nfig + 1;
figure(nfig)
LL = (Kx2.^2 + Kz2.^2).^(1/2);
hs=quiver(X1, Z, Kx2./LL, Kz2./LL, 0.5);  
set(hs,'AutoScale', 'on','AutoScaleFactor', 2);
ylabel('y axis')
xlabel('x axis')
title('Surface Current K(x,y) with quiver')
axis equal square
grid on
hold off

nfig = nfig + 1;
figure(nfig)
hs=streamslice(X1, Z, Kx2, Kz2, 2); 
set(hs,'Color','m','Linewidth',0.8)
ylabel('y axis')
xlabel('x axis')
title('Surface Current K(x,y)with streamslice')
axis equal
grid on
hold off
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%                                                  %
%                 Functions                        %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

function A_x = Cal_Of_Vec_Pot_x(x,y,z,d,h, mu_0, m0)
  R1 = ((x - d).^2 + (y - h).^2 + z.^2).^(3/2);
  R2 = ((x - d).^2 + (y + h).^2 + z.^2).^(3/2);
  R3 = ((x + d).^2 + (y - h).^2 + z.^2).^(3/2);
  R4 = ((x + d).^2 + (y + h).^2 + z.^2).^(3/2);

  factor = (mu_0 * m0)/(4*pi);
  A_x = factor.*(z./R1 - z./R2 + z./R3 - z./R4);
end
function A_z = Calc_of_Vec_Pot_z(x,y,z,d,h, mu_0, m0)
  R1 = ((x - d).^2 + (y - h).^2 + z.^2).^(3/2);
  R2 = ((x - d).^2 + (y + h).^2 + z.^2).^(3/2);
  R3 = ((x + d).^2 + (y - h).^2 + z.^2).^(3/2);
  R4 = ((x + d).^2 + (y + h).^2 + z.^2).^(3/2);
  factor = (mu_0 * m0)/(4*pi);
  A_z = factor.*(-(x-d)./R1 + (x-d)./R2 - (x+d)./ R3 + (x+d)./R4);
end

%magnetic-field
function H_x = Calc_Magnetic_Field_x(x,y,d,h,m0)
  R1 = ((x - d).^2 + (y - h).^2) .^(1/2);
  R2 = ((x - d).^2 + (y + h).^2).^(1/2);
  R3 = ((x + d).^2 + (y - h).^2).^(1/2);
  R4 = ((x + d).^2 + (y + h).^2).^(1/2);
  factor = (m0)/(4*pi);
  H_x = factor * ( ((3*(y-h).*(x-d))./(R1.^5)) - ((3*(y+h).*(x-d))./(R2.^5)) + ((3*(y-h).*(x+d))./(R3.^5)) - ((3*(y+h).*(x+d))./(R4.^5)));
end
function H_y = Calc_Magnetic_Field_y(x,y,d,h,m0)
  R1 = ((x - d).^2 + (y - h).^2) .^(1/2);
  R2 = ((x - d).^2 + (y + h).^2).^(1/2);
  R3 = ((x + d).^2 + (y - h).^2).^(1/2);
  R4 = ((x + d).^2 + (y + h).^2).^(1/2);
  factor = (m0)/(4*pi);
  H_y = factor * ( (((3*(y-h).^2) ./(R1.^5))- 1./(R1.^3)) + ((((-3)*(y+h).^2)./(R2.^5)) + 1./(R2.^3)) + (((3*(y-h).^2)./(R3.^5)) - 1./(R3.^3)) + ((((-3)*(y+h).^2)./(R4.^5)) + 1./(R4.^3)));
end

%surface-current
function k_y = Calc_Surface_Current_y(x,y,z,d,h,m0)
  R13 = ((x - d).^2 + (y - h).^2 +z.^2) .^(1/2);
  R24 = ((x - d).^2 + (y + h).^2 +z.^2).^(1/2);
  factor = (3*m0.*z)./(2*pi);
  k_y = factor * ((-y+h)./(R13.^5) + (y+h)./(R24.^5));
end
function k_z = Calc_Surface_Current_z(x,y,z,d,h,m0)
  R13 = ((x - d).^2 + (y - h).^2) .^(1/2);
  R24 = ((x - d).^2 + (y + h).^2).^(1/2);
  factor = m0/(2*pi);
  k_z = factor * ( (3*(y-h).^2 ./(R13.^5) - 1./(R13.^3)) +(-3*(y+h).^2 ./(R24.^5) + 1./(R24.^3)) );
end

function k2_x = Calc_Surface_Current_new_x(x,y,z,d,h,m0)
    R12 = ((x - d).^2 + (y - h).^2 + z.^2).^(1/2);          
    R34 = ((x + d).^2 + (y - h).^2 + z.^2).^(1/2);
    factor = ((-3*m0*h)/(2*pi));
    k2_x = factor * (z./(R34.^5) - z./(R12.^5) );
end
function k2_z = Calc_Surface_Current_new_z(x,y,z,d,h,m0)
    R12 = ((x - d).^2 + (y - h).^2 + z.^2).^(1/2);
    R34 = ((x + d).^2 + (y - h).^2 + z.^2).^(1/2);
    factor = 3*m0*h/(2*pi);
    k2_z = factor * (((x+d)./(R34.^5)) - ((x-d)./(R12.^5)) );
end
%end