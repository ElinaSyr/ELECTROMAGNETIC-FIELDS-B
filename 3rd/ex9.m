%data
a=0.1;
I1 = 1;
I2 = 1;
nfig = 0;
m0 = 4*pi*1e-7;

nfig = question_d(a,a,I1,I2,nfig); %ερωτημα δ για α = h
%nfig = question_d(a/2,a,I1,I2,nfig);ερωτημα δ για α/2 = h
%nfig = question_d(a/4,a,I1,I2,nfig); ερωτημα δ για α/4 = h

nfig = question_e(a/2,a,I1,I2,nfig);
nfig = question_z(0.025,nfig);
%nfig = question_z(0.05,nfig); % ζ ερώτημα για h = 0.05 
%nfig = question_z(0.1,nfig); %ζ ερώτημα για h = 0.1

nfig = question_st(h,a,I,nfig);

function ret = question_d(h,a,I1,I2,nfig)
    %grid
    Npoints = 250;
    zmin = -2*h;
    zmax = 2*h;
    zz = zmin:(zmax-zmin)/Npoints:zmax;

    %
    Hz = magneticField(zz, I1, I2, a, h);

    nfig = nfig+1;
    hold off;
    figure(nfig)
    plot(zz, Hz)
    
    %differentiation
    diff_Hz = gradient(Hz);
    for i=1:4
        nfig = nfig+1;
        hold off;
        figure(nfig)
        plot(zz, diff_Hz)
        diff_Hz = gradient(diff_Hz);
    end
   
    %functions
    function H_z = magneticField(z, I1, I2, a, h) 
        H_z = (I1 * a^2) ./ (2*(a.^2+(z-h).^2).^(3/2)) + (I2 * a^2) ./ (2*(a.^2+(z+h).^2).^(3/2));
    end

    ret = nfig;
end



function ret = question_e(h,a,I1,I2,nfig)
    %grid
    xmin = -3*a;
    xmax = 3*a;
    zmin = -(h+2*a);
    zmax = h+2*a;
    
    Npoints =251;
    x1 = xmin: (xmax-xmin)/Npoints:xmax;
    z1 = zmin: (zmax-zmin)/Npoints:zmax;
    
    [X,Z] = meshgrid(x1,z1);
    A_array = zeros(length(x1),length(z1));
    %
   
    for ix=1:length(x1)
        for iz=1:length(z1)
            xx1 = X(ix, iz);
            zz1 = Z(ix, iz);
            [func] = integral_Of_A(xx1,0,zz1,a,h,0,2*pi);
            A_array(ix, iz) = (I1*a/(4*pi))*func;
        end
     end
    
    %surface
    nfig = nfig+1;
    figure(nfig);
    hold off
    surface(X,Z,A_array), shading interp
    axis equal
    colorbar
    
    nfig = nfig + 1;
    figure(nfig)
    hold off
    cont = [0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.175, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,...
    1.25, 1.5, 2.0, 3.0, 5.0];
    [CS,H] = contour(X, Z, A_array,cont,'Linewidth',1,'Color','b');
    clabel(CS, H, cont);
    xlabel('x')
    ylabel('z')
    title('Vector Potential with contour function')

    function A = integral_Of_A(x,y,z,a,h,umin,umax)
            
            A1 = integral(@(u) (cos(u))/((x.^2 + y.^2 + (z-h).^2 + a.^2 -2*a.*((x.^2 + y.^2)^(1/2))*cos(u)).^(1/2)), umin, umax, 'ArrayValued', true);
            A2 = integral(@(u) (cos(u))/((x.^2 + y.^2 + (z+h).^2 + a^2 -2*a.*((x.^2 + y.^2).^(1/2))*cos(u)).^(1/2)), umin, umax, 'ArrayValued', true);

            A = A1 + A2;
    end
    ret = nfig;
end
function ret = question_st(h,a,I,nfig)
    xmin = -3*a;
    xmax = 3*a;
    zmin = -(h+2*a);
    zmax = h+2*a;
    
    Npoints =251;
    x1 = xmin: (xmax-xmin)/Npoints:xmax;
    z1 = zmin: (zmax-zmin)/Npoints:zmax;
    
    [X,Z] = meshgrid(x1,z1);
    
     Mag_Field_x = zeros (length(x1),length(z1));
     Mag_Field_z = zeros (length(x1),length(z1));
     
     for ix=1:length(x1)
        for iz=1:length(z1)
            xx1 = X(ix, iz);
            zz1 = Z(ix, iz);
            [func1, func2] = Calc_Of_Magn_Field(xx1,0,zz1,a,h,I);
            Mag_Field_x(ix, iz) = func1;
            Mag_Field_z(ix, iz) = func2;
        end
     end
     
     
    nfig = nfig + 1;
    figure(nfig)
    LL = (Mag_Field_x.^2 + Mag_Field_z.^2).^(1/2);
    quiver(X, Z, Mag_Field_x./LL, Mag_Field_z./LL, 0.5);
    ylabel('z')
    xlabel('x')
    title('Magnetic Field H(x,z) with quiver')
    axis equal
    grid on
    hold off

    nfig = nfig + 1;
    figure(nfig)
    hs=streamslice(X, Z, Mag_Field_x, Mag_Field_z, 2); 
    set(hs,'Color','m','Linewidth',0.8)
    ylabel('z')
    xlabel('x')
    title('Magnetic Field H(x,z) with streamslice')
    axis equal
    grid on
    hold off
%Με υπολογισμούς του β ερωτήματος προκύπτουν τα παρακάτω 
   function [Hx, Hz] = Calc_Of_Magn_Field(x,y,z,a,h,I)
    
        RR1 = (x^2 + y^2 +(z-h)^2)^(1/2);
        A1_con = RR1^2 + a^2 - 2*a*((x^2 + y^2)^(1/2));
        B1_con = RR1^2 + a^2 + 2*a*((x^2 + y^2)^(1/2));
        LL1 = 1-(A1_con/B1_con);
        [K1, E1] = ellipke(LL1);
        
        RR2 = (x^2 + y^2 +(z+h)^2)^(1/2);
        A2_con = RR2^2 + a^2 - 2*a*((x^2 + y^2)^(1/2));
        B2_con = RR2^2 + a^2 + 2*a*((x^2 + y^2)^(1/2));
        LL2 = 1-(A2_con/B2_con);
        [K2, E2] = ellipke(LL2);
        
        Hx1 = ((I/pi)*a^2*(z-h)*(1/RR1)*(1/(A1_con*((B1_con)^(1/2))))*E1)*(x/RR1) ...
            + ((I/pi)*a^2*(z+h)*(1/RR2)*(1/(A2_con*((B2_con)^(1/2))))*E2)*(x/RR2);
        Hz1 = ((I/pi)*a^2*(z-h)*(1/RR1)*(1/(A1_con*((B1_con)^(1/2))))*E1)*((z-h)/RR1)...
            + ((I/pi)*a^2*(z+h)*(1/RR2)*(1/(A2_con*((B2_con)^(1/2))))*E2)*((z+h)/RR2);


        Hx2 = ((I/pi)*RR1*((RR1^2  + (a^2/RR1^2)*((z-h)^2 - x^2 - y^2))*E1-A1_con*K1)/(2*A1_con*((B1_con)^(1/2))*((x^2 + y^2)^(1/2))))*((z-h)*x/(RR1*((x^2+y^2)^(1/2)))) +...
               ((I/pi)*RR2*((RR2^2  + (a^2/RR2^2)*((z+h)^2 - x^2 - y^2))*E2-A2_con*K2)/(2*A2_con*((B2_con)^(1/2))*((x^2 + y^2)^(1/2))))*((z+h)*x/(RR2*((x^2+y^2)^(1/2))));


        Hz2 = ((I/pi)*RR1*((RR1^2  + (a^2/RR1^2)*((z-h)^2 - x^2 - y^2))*E1-A1_con*K1)/(2*A1_con*((B1_con)^(1/2))*((x^2 + y^2)^(1/2))))*(-((x^2+y^2)^(1/2))/RR1) +...
               ((I/pi)*RR2*((RR2^2  + (a^2/RR2^2)*((z+h)^2 - x^2 - y^2))*E2-A2_con*K2)/(2*A2_con*((B2_con)^(1/2))*((x^2 + y^2)^(1/2))))*(-((x^2+y^2)^(1/2))/RR2);
        
        Hx = Hx1 + Hx2;
        Hz = Hz1 + Hz2;
    end
    ret = nfig; 
end 

function ret = question_z(h,nfig)
    amin = 0;
    amax = 0.25; 
    Npoints = 200;
    a = amin:(amax - amin)/Npoints:amax; 
    k = 2.*a.*(1./(((2*h).^2 + 4*a.^2).^(1/2)));
    [K, E] = ellipke(k);
    L_coil = a.*((2-k.^2).*K./k - 2*E./k);
    
  
       
   nfig = nfig + 1;
   figure(nfig)
   plot(a, L_coil)
   hold on
   xlabel("Radious a")
   ylabel("Induction ")
   title("Induction of the coils for different radious ")
   

   ret = nfig;
end