function Cmat = Cmat_transvIsoV3b(ndir, eps, E_m, nu_m, Ax, Bx, theta0, E_f, nu_f, v_f)
    %Note: A and B variable names changed to Ax Bx 
    e11 = eps(1);
    e22 = eps(2);
    e12 = eps(3);
    
    e33 = 0;
    e23 = 0;
    e13 = 0;
        
 if ndir == 1
    C = zeros(6);
    e11 = eps(1);
    e22 = eps(2);
    e12 = eps(3);
    C(1,1) = - v_f*((4*Ax*sin(theta0/2)^2)/(sin(theta0/2)^2*(e11 + 1)^2 - 1) - (Bx*((theta0*sin(theta0/2)^2)/(asin(sin(theta0/2)*(e11 + 1))^2*(sin(theta0/2)^2*(e11 + 1)^2 - 1)) - (sin(theta0/2)^3*(2*e11 + 2))/(1 - sin(theta0/2)^2*(e11 + 1)^2)^(3/2) + (theta0*sin(theta0/2)^3*(2*e11 + 2))/(2*asin(sin(theta0/2)*(e11 + 1))*(1 - sin(theta0/2)^2*(e11 + 1)^2)^(3/2))))/theta0 + (Ax*sin(theta0/2)^3*(2*e11 + 2)*(theta0 - 2*asin(sin(theta0/2)*(e11 + 1))))/(1 - sin(theta0/2)^2*(e11 + 1)^2)^(3/2) + (2*E_f*e22^2*sin(theta0/2)^2*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^4*(E_m/E_f - 1)*(nu_f - 1))/((sin(theta0/2)^2*(e11 + 1)^2 - 1)*(2*nu_f^2 + nu_f - 1)) - (6*E_f*e22^2*sin(theta0/2)^2*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^2*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^2*(E_m/E_f - 1)*(nu_f - 1))/((sin(theta0/2)^2*(e11 + 1)^2 - 1)*(2*nu_f^2 + nu_f - 1)) + (E_f*e22^2*sin(theta0/2)^3*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^3*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*(2*e11 + 2)*(E_m/E_f - 1)*(nu_f - 1))/((1 - sin(theta0/2)^2*(e11 + 1)^2)^(3/2)*(2*nu_f^2 + nu_f - 1)) - (2*E_f*e22*e33*nu_f*sin(theta0/2)^2*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^2*(E_m/E_f - 1))/((2*nu_f - 1)*(sin(theta0/2)^2*(e11 + 1)^2 - 1)*(nu_f + 1)) + (2*E_f*e22*e33*nu_f*sin(theta0/2)^2*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^2*(E_m/E_f - 1))/((2*nu_f - 1)*(sin(theta0/2)^2*(e11 + 1)^2 - 1)*(nu_f + 1)) - (E_f*e22*e33*nu_f*sin(theta0/2)^3*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*(2*e11 + 2)*(E_m/E_f - 1))/((2*nu_f - 1)*(1 - sin(theta0/2)^2*(e11 + 1)^2)^(3/2)*(nu_f + 1))) - ((2*E_m)/(2*nu_m + 2) - (E_m*nu_m)/((2*nu_m - 1)*(nu_m + 1)))*(v_f - 1);
    C(1,2) = (E_m*nu_m*(v_f - 1))/((2*nu_m - 1)*(nu_m + 1)) - v_f*((4*E_f*e22*sin(theta0/2)*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^3*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*(E_m/E_f - 1)*(nu_f - 1))/((1 - sin(theta0/2)^2*(e11 + 1)^2)^(1/2)*(2*nu_f^2 + nu_f - 1)) - (2*E_f*e33*nu_f*sin(theta0/2)*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*(E_m/E_f - 1))/((2*nu_f - 1)*(1 - sin(theta0/2)^2*(e11 + 1)^2)^(1/2)*(nu_f + 1)));
    C(2,2) = -1/((v_f - 1)/((2*E_m)/(2*nu_m + 2) - (E_m*nu_m)/((2*nu_m - 1)*(nu_m + 1))) + (v_f*(2*nu_f^2 + nu_f - 1))/(E_f*(cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^4*(E_m/E_f - 1) - E_m/E_f)*(nu_f - 1)));
    C(2,1) = C(1,2);
    C(6,6) = -1/(((2*nu_m + 2)*(v_f - 1))/(2*E_m) - (v_f*(2*nu_f + 2))/E_f);
    
    Cmat = [C(1,1) C(1,2) 0;...
         C(2,1) C(2,2) 0;...
         0 0 C(6,6)];
 end
 
 if ndir == 2   
    e11 = eps(2);
    e22 = eps(1);
    e12 = eps(3);
    C = zeros(6);   
    C(2,2) = - v_f*((4*Ax*sin(theta0/2)^2)/(sin(theta0/2)^2*(e11 + 1)^2 - 1) - (Bx*((theta0*sin(theta0/2)^2)/(asin(sin(theta0/2)*(e11 + 1))^2*(sin(theta0/2)^2*(e11 + 1)^2 - 1)) - (sin(theta0/2)^3*(2*e11 + 2))/(1 - sin(theta0/2)^2*(e11 + 1)^2)^(3/2) + (theta0*sin(theta0/2)^3*(2*e11 + 2))/(2*asin(sin(theta0/2)*(e11 + 1))*(1 - sin(theta0/2)^2*(e11 + 1)^2)^(3/2))))/theta0 + (Ax*sin(theta0/2)^3*(2*e11 + 2)*(theta0 - 2*asin(sin(theta0/2)*(e11 + 1))))/(1 - sin(theta0/2)^2*(e11 + 1)^2)^(3/2) + (2*E_f*e22^2*sin(theta0/2)^2*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^4*(E_m/E_f - 1)*(nu_f - 1))/((sin(theta0/2)^2*(e11 + 1)^2 - 1)*(2*nu_f^2 + nu_f - 1)) - (6*E_f*e22^2*sin(theta0/2)^2*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^2*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^2*(E_m/E_f - 1)*(nu_f - 1))/((sin(theta0/2)^2*(e11 + 1)^2 - 1)*(2*nu_f^2 + nu_f - 1)) + (E_f*e22^2*sin(theta0/2)^3*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^3*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*(2*e11 + 2)*(E_m/E_f - 1)*(nu_f - 1))/((1 - sin(theta0/2)^2*(e11 + 1)^2)^(3/2)*(2*nu_f^2 + nu_f - 1)) - (2*E_f*e22*e33*nu_f*sin(theta0/2)^2*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^2*(E_m/E_f - 1))/((2*nu_f - 1)*(sin(theta0/2)^2*(e11 + 1)^2 - 1)*(nu_f + 1)) + (2*E_f*e22*e33*nu_f*sin(theta0/2)^2*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^2*(E_m/E_f - 1))/((2*nu_f - 1)*(sin(theta0/2)^2*(e11 + 1)^2 - 1)*(nu_f + 1)) - (E_f*e22*e33*nu_f*sin(theta0/2)^3*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*(2*e11 + 2)*(E_m/E_f - 1))/((2*nu_f - 1)*(1 - sin(theta0/2)^2*(e11 + 1)^2)^(3/2)*(nu_f + 1))) - ((2*E_m)/(2*nu_m + 2) - (E_m*nu_m)/((2*nu_m - 1)*(nu_m + 1)))*(v_f - 1);
    C(2,1) = (E_m*nu_m*(v_f - 1))/((2*nu_m - 1)*(nu_m + 1)) - v_f*((4*E_f*e22*sin(theta0/2)*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^3*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*(E_m/E_f - 1)*(nu_f - 1))/((1 - sin(theta0/2)^2*(e11 + 1)^2)^(1/2)*(2*nu_f^2 + nu_f - 1)) - (2*E_f*e33*nu_f*sin(theta0/2)*cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*sin(pi/2 - asin(sin(theta0/2)*(e11 + 1)))*(E_m/E_f - 1))/((2*nu_f - 1)*(1 - sin(theta0/2)^2*(e11 + 1)^2)^(1/2)*(nu_f + 1)));
    C(1,1) = -1/((v_f - 1)/((2*E_m)/(2*nu_m + 2) - (E_m*nu_m)/((2*nu_m - 1)*(nu_m + 1))) + (v_f*(2*nu_f^2 + nu_f - 1))/(E_f*(cos(pi/2 - asin(sin(theta0/2)*(e11 + 1)))^4*(E_m/E_f - 1) - E_m/E_f)*(nu_f - 1)));
    C(1,2) = C(2,1);
    C(6,6) = -1/(((2*nu_m + 2)*(v_f - 1))/(2*E_m) - (v_f*(2*nu_f + 2))/E_f);
    
    Cmat = [C(1,1) C(1,2) 0;...
         C(2,1) C(2,2) 0;...
         0 0 C(6,6)];
 end
% % Test linear elastic material 
% E=1e5; v=1/3;
% Cmat = (E/((1+v)*(1-2*v)))*[1-v v 0;v 1-v 0; 0 0 (1-2*v)/2];