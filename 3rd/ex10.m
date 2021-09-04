a=1; 
M0=1;
m0=4*pi*1e-7;
nfig = 0;

%Grid 
xmin = -3*a;
xmax = 3*a; 
ymin = -3*a;
ymax = 3*a; 
Npoints = 151; 
x1 = xmin: (xmax -xmin)/Npoints : xmax;
y1 = ymin: (ymax -ymin)/Npoints : ymax;

[X,Y]=meshgrid(x1,y1);

A_array = zeros(length(x1),length(y1)); 

%Vector Potential 
for ix=1:length(x1)
    for iy=1:length(y1)
        xx1=X(ix,iy);
        yy1=Y(ix,iy);
        [theta,rho,Z]=cart2pol(xx1,yy1,1); %μετατροπή από καρτεσιανές σε κυλινδρικές 
        A_array(ix,iy)=-a*integral(@(phi2) vector_pot(phi2,rho,a,theta),0,2*pi);
    end
end


%Surface + shading interp 
nfig = nfig + 1; 
figure(nfig)
hold off
surface(X,Y,A_array),shading interp
hold on
axis equal
title('Vector Potential with surface function');
xlabel('x axis')
ylabel('y axis')

% Contour
nfig = nfig + 1;
figure(nfig)
cont=2.7*(-0.9:0.1:0.9);
hold off
[CS,H]=contour(X,Y,A_array,cont,'LineWidth',1);
clabel(CS,H,cont);
axis equal
hold on
title('Vector Potential with contour function')
xlabel('x axis')
ylabel('y axis')

%αρχικοποιήση 
Bx = zeros (length(x1),length(y1));
By = zeros (length(x1),length(y1));
Hx = zeros (length(x1),length(y1));
Hy = zeros (length(x1),length(y1));
% B AND H 
for ix=1:length(x1)
    for iy=1:length(y1)
        xx1=X(ix,iy);
        yy1=Y(ix,iy);
        [theta,rho,Z]=cart2pol(xx1,yy1,1); 
        factor =(m0*M0/(2*pi));
        
        Bx(ix,iy)= factor*a*integral(@(phi2) mag_ind_x(phi2,rho,a,theta,yy1),0,2*pi);
        By(ix,iy)=-factor*a*integral(@(phi2) mag_ind_y(phi2,rho,a,theta,xx1),0,2*pi);
        
        if (rho>a)
            Hx(ix,iy)=Bx(ix,iy)/m0;
            Hy(ix,iy)=By(ix,iy)/m0;
        else
            Hx(ix,iy)=Bx(ix,iy)/m0;
            Hy(ix,iy)=By(ix,iy)/m0-M0;
        end
    end
end


nfig = nfig + 1; 
figure(nfig)
hold off 
LL= (Bx.^2+By.^2).^(1/2);
quiver(X,Y,Bx./LL,By./LL,0.5);
axis equal
hold on

nfig = nfig + 1; 
figure(nfig)
hs=streamslice(X,Y,Bx,By,1);
title('Magnetic Induction B(x,y)')
xlabel('x axis')
ylabel('y axis')



nfig = nfig + 1; 
figure(nfig)
hold off
LL=(Hx.^2+Hy.^2).^(1/2);
quiver(X,Y,Hx./LL,Hy./LL,0.5);

axis equal
hold on

nfig = nfig + 1; 
figure(nfig)
hs=streamslice(X,Y,Hx,Hy,1);
title('Magnetic field H(x,y)')
xlabel('x axis')
ylabel('y axis')
hold off

function AA = vector_pot(phi2,r,a,phi1)
AA=cos(phi2).*log(a./(r.^2 +a^2 -2.*r.*a*cos(phi1-phi2)).^(1/2));
end
function B_x = mag_ind_x(phi2,a,rt,phi1,yy)
B_x=((yy-a.*sin(phi2)).*cos(phi2))./(rt.^2+a^2-2.*rt.*a.*cos(phi1-phi2));
end
function B_y = mag_ind_y(phi2,a,rt,phi1,xx)
B_y=((xx-a.*cos(phi2)).*cos(phi2))./(rt.^2+a^2-2.*rt.*a.*cos(phi1-phi2));
end