clc
clear all 
close all 

%% Parameters of Blasius Equation
U_inf = 1;
L = 10;
mu = 1.789E-5;
rho = 1.225;
nu = mu/rho;
A = sqrt(nu/U_inf);
h = 0.01;
tol = 1.0e-5;

%% Numerical Solution of Blasius Equation Using Runge-Kutta

zeta_old2 = 0.2;
zeta_old1 = 0.3;

[f1,y1,z1] = Runge_Kutta4(zeta_old2);

[f2,y2,z2] = Runge_Kutta4(zeta_old1);

n = length(z1);

yb_old2 = y1(n)
yb_old1 = y2(n)
beta = 1;

while (abs(yb_old1 - beta) >= tol)
    
    zeta = zeta_old1 - (yb_old1 - beta)*(zeta_old1 - zeta_old2)/(yb_old1 - yb_old2)
    aux = zeta_old1
    zeta_old1 = zeta
    zeta_old2 = aux
    
    yb_old2 = yb_old1
    [f2,y2,z2] = Runge_Kutta4(zeta_old1);
    yb_old1 = y2(n)

end
