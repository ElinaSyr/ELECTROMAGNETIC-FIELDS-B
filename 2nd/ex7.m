 L = 0.5;
 h = 0.5;
 e0 = 8.854e-12;
 nfig = 1;
 a = 0.0025; 
 N = 50; 
 Dx = 2*L/N;
 range = linspace(-4,4,200);
 [x,z] = meshgrid(range,range);
 factor1 = 1/(4*pi); 
 factor2 = -2/3; 
 factor3 = 1/(12*pi);
 %Area1 = z>h;
 %Area2 = (z<h).*(z>=0);
 Area3 = z<0;
 Area4 = z>0;
 Phi1a = factor1.*log(((L-x)+sqrt((L-x).^2 + (z-h).^2))./(-L-x + sqrt((L+x).^2 + (z-h).^2)));
 Phi2a = factor1.*factor2.*log((L-x + sqrt((L-x).^2 + (z+h).^2))./(-L-x + sqrt((L+x).^2 + (z+h).^2)));
 Phi3b = factor3.*log((L-x + sqrt((L-x).^2 + (z-h).^2))./(-L-x + sqrt((L+x).^2 + (z-h).^2)));
 Phi1=Phi1a+Phi2a;
 Phi = Phi1.*Area4+Phi3b.*Area3;
 figure(nfig);
 colorbar; hold on 
 surface(x,z,Phi,'LineStyle','None'), shading interp; 
 title('Potential (x,z)/ 0, h = 1m, 1/ 0 = 1, 2/ 0 = � � � � � � 5','Fontsize',12,'FontWeight','bold')
 xlabel('x-codinate','Fontsize',12,'FontWeight','bold')
 ylabel('y-cordinate','Fontsize',12,'FontWeight','bold')
 cont =[0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.2 0.3 0.4 0.5];
 nfig = nfig+1;
 figure(nfig);
 hold off
 [CS,H] = contour(x,z,Phi,cont,'Linewidth',1);
 clabel(CS,H,cont);
 axis equal
 title('Equipotential Lines, h = 1m, e1/e0 = 1, e2/e0 = 5','Fontsize',12,'FontWeight','bold')
 xlabel('x-cordinate','Fontsize',12,'FontWeight','bold')
 ylabel('z-cordinate','Fontsize',12,'FontWeight','bold')
Dx=4/200;
Dz=4/200;
[Exarray1, Eyarray1] = gradient(Phi,Dx,Dz);
Exarray = - Exarray1;
Eyarray = - Eyarray1;
%Plot Electric Field 
nfig = nfig+1;
figure(nfig)
hold off
LL = ((Exarray).^2 + (Eyarray).^2).^(1/2);
axis([-2 2 -2 2]); hold on;
quiver(x,z,Exarray./LL,Eyarray./LL,0.5); 
title ('Electrical Field','Fontsize',12,'FontWeight','bold')
xlabel('x-cordinate','Fontsize',12,'FontWeight','bold')
ylabel('z-cordinate','Fontsize',12,'FontWeight','bold')
%%% Some other constants %%%
 N = 50; 
 Dx = 2*L/N;
 dx =4/100;
 range_n = -2:dx:2;
 x = -L:Dx:L; 
 %% Factors %%
 e1 = e0 ;
 e2 = 5 * e0; 
 factor = 1/(4*pi*e1) ;
 factor_1 = -2/3 ; %%%d
for w=1:N
 xx(w)= (1/2)*(x(w+1)+x(w));
end
for i = 1:N
 for j = 1:N
 Rij(d,l) = abs(xx(d) - xx(l));
 Ri_j(d,l) = ((xx(d) - xx(l)).^2 + 4*h^2).^(1/2);
 if i==j
 Aij(i,j) = factor*(2*log(Dx/a)+ factor_1 * Dx/(2*h));
 else
 Aij(i,j) = factor*Dx.* (1./Rij(i,j) + factor_1 .* 1./Ri_j(i,j));
 end
 end
end 
ones_vector = ones(N); %% created a vector of 1 column and N lines 
Phi0 = 1; 
lambda0= (Aij\ones_vector) .*Phi0;
lambda = (1./(Dx.*Sum(lambda0(:,1)))) .* (Aij\ones_vector) .*Phi0;
l_n = lambda';
nfig = nfig+1;
figure(nfig);
hold off
plot(x,l_n);
xlabel('Normalized Length of Conductive wire','Fontsize',12,'FontWeight','bold')
 ylabel('Normalized Line Charge Density, (x)/0 (C/m) � ','Fontsize',12,'FontWeight','bold')
axis equal;
hold off
%Potential on the surface
factor_4 = Dx/(2*pi*(e1 + e2));
for ix = 1:length(x)
 Sum=0;
 for j = 1:N
 R = ((x(j) - xx(ix)).^2 + h^2).^(1/2);
 Sum = Sum + factor_4 * lambda(j)/R;
 end
 Volt_a(ix) = Sum;
end
nfig = nfig+1;
figure(nfig);
hold off
plot(xx,Volt_a*e0);