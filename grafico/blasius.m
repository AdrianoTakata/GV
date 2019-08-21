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

%% Numerical Solution of Blasius Equation Using Runge-Kutta
f1 = @(x, y1, y2, y3) y2;
f2 = @(x, y1, y2, y3) y3;
f3 = @(x, y1, y2, y3) -y1*y3;
eta = 0:h:10;
x = 0:h:10;
y1(1) = 0;
y2(1) = 0;
y3(1) = 0.4696;
for i = 1:(length(eta)-1)
a = h.*[f1(eta(i), y1(i), y2(i), y3(i)), f2(eta(i), y1(i), y2(i), y3(i)), f3(eta(i), y1(i), y2(i), y3(i))];
b = h.*[f1(eta(i), y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2), f2(eta(i)+h/2, y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2), f3(eta(i)+h/2, y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2)];
c = h.*[f1(eta(i), y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2), f2(eta(i)+h/2, y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2), f3(eta(i)+h/2, y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2)];
d = h.*[f1(eta(i), y1(i)+c(1), y2(i)+c(2), y3(i)+c(3)), f2(eta(i)+h, y1(i)+c(1), y2(i)+c(2), y3(i)+c(3)), f3(eta(i)+h, y1(i)+c(1), y2(i)+c(2), y3(i)+c(3))];
y3(i+1) = y3(i)+ 1/6*(a(3)+2*b(3)+2*c(3)+d(3));
y2(i+1) = y2(i)+ 1/6*(a(2)+2*b(2)+2*c(2)+d(2));
y1(i+1) = y1(i)+ 1/6*(a(1)+2*b(1)+2*c(1)+d(1));
end

%% Plotting and Visualization
figure(1)
plot(y1, eta, y2, eta, y3, eta, 'LineWidth', 2)
xlim([0 2])
title('Solution of Blasius eqution', 'FontSize', 14);
xlabel('f, f'' and f''''', 'FontSize', 20);
ylabel('\eta', 'FontSize', 20);
grid on
Legend1 = {'f(\eta)', 'f''(\eta)', 'f''''(\eta)'};
legend(Legend1, 'FontSize', 14);

%% Velocity Profile and the BL thickness distribution
figure(2)
hold all
for i = 1:length(x)
delta(i) = 5*sqrt(x(i))*A;
end
plot(x, delta, 'LineWidth', 2, 'color', 'black');


position = [1 4 5];
for j = 1:length(position)
y{j} = eta*sqrt(position(j))*A;
plot(y2+position(j), y{j}, 'LineWidth', 2)
end
Legend2 = {'\delta (x)', 'Velocity Profile at x = 1', 'Velocity Profile at x = 4', 'Velocity Profile at x = 5'};
title('Velocity profiles at x = 1, 4 and 5 plus BL thickness distribution along the plate', 'FontSize', 14);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
legend(Legend2, 'FontSize', 14);
ylim([0, 2*max(y{j})])
grid on

%% shear stress at the wall as function of x
figure(3)
for i = 1 : length(x)
tau_x(i) = mu*U_inf*y3(1)*sqrt(U_inf/2/nu/x(i));
end
plot(x, tau_x, 'LineWidth', 2, 'color', 'black')
title(' shear stress at the wall along the flat plate', 'FontSize', 14);
xlabel('x', 'FontSize', 20);
ylabel('\tau_w_a_l_l', 'FontSize', 20);
axis tight
grid on

%% shear stress as function of eta at x = 1, 4, 5
figure(4)
hold all
for i = 1:length(position)
tau_eta = mu*U_inf*y3*sqrt(U_inf/2/nu/position(i));
plot(x, tau_eta, 'LineWidth', 2)
end
title(' shear stress as funcion of \eta ', 'FontSize', 14);
xlabel('\eta', 'FontSize', 20);
ylabel('\tau', 'FontSize', 20);
axis tight
grid on
Legend3 = {'\tau(\eta) at x = 1', '\tau(\eta) at x = 4', '\tau(\eta) at x = 5'};
legend(Legend3, 'FontSize', 14);


%% shear derivative as function of eta at x = 1, 4, 5
figure(5)
hold all
y4 = -1/2*y1.*y3;
for i = 1 : length(position)
tau_d = mu*U_inf^2/2/nu/position(i)*y4; %shear derivative partial(tau)/partial(y) at the wall y=0
plot(eta, tau_d, 'LineWidth', 2)
end
title(' shear stress derivative as function of \eta ', 'FontSize', 14);
xlabel('\eta', 'FontSize', 20);
ylabel('\tau''', 'FontSize', 20);
axis tight
grid on
Legend4 = {'\tau''(\eta) at x = 1', '\tau''(\eta) at x = 4', '\tau''(\eta) at x = 5'};
legend(Legend4, 'FontSize', 14);


%% local shear cofficient as a function of x
figure(6)
hold all
for i = 1: length(x)
cfx(i) = tau_x(i)*2/rho/U_inf^2;
% cfxt(i) = 0.664/sqrt(rho*U_inf*x(i)/mu);
end
plot(x, cfx, 'LineWidth', 2, 'color', 'black')
title(' Local shear coefficient along the flat plate', 'FontSize', 14);
xlabel('x', 'FontSize', 20);
ylabel('C_f_x', 'FontSize', 20);
axis tight
grid on
%% Total Skin friction coefficient
CF_T = 0;
dx = x(2)-x(1);
for i = 2 : length(x)
CF_T = CF_T + 1/L*cfx(i)*dx;
end
CDF_T = CF_T*2;
%% Accumulative Drag Force as function of x
CF = 0;
for i = 2 : length(x)
CF(i) = 1/L*cfx(i)*dx;
CF(i) = CF(i-1)+CF(i);
CDF(i) = CF(i)*2;
D(i) = 1/2*rho*U_inf^2*L*CDF(i);
end

figure(7)
plot(x, CF, 'LineWidth', 2, 'color', 'black');
title(' Accumulative Total friction coeff', 'FontSize', 14);
xlabel('x', 'FontSize', 20);
ylabel('C_F', 'FontSize', 20);
axis tight
grid on

figure(8)
plot(x, D, 'LineWidth', 2, 'color', 'black');
title(' Accumulative Drag Force along the flat plate', 'FontSize', 14);
xlabel('x', 'FontSize', 20);
ylabel('Drag Force (D)', 'FontSize', 20);
axis tight
grid on

figure(9)
y1 = zeros(1, length(x));
position1 = ones(1, length(x))*position(1);
position2 = ones(1, length(x))*position(2);
position3 = ones(1, length(x))*position(3);
C1 = ones(1, (length(x)+1))*position(1);
C2 = ones(1, (length(x)+1))*position(2);
C3 = ones(1, (length(x)+1))*position(3);

hold all
fill3([position1, position1(1)], [y2, 0], [y{1}, y{1}(length(x))], C1, 'LineWidth', 2)
fill3([position2, position2(2)], [y2, 0], [y{2}, y{2}(length(x))], C2, 'LineWidth', 2)
fill3([position3, position3(3)], [y2, 0], [y{3}, y{3}(length(x))], C3, 'LineWidth', 2)
plot3(x, y1, delta, 'LineWidth', 2, 'color', 'black')
grid on
title('BL thickness and velocity profile of the flow', 'FontSize', 14);
xlabel('x', 'FontSize', 20);
ylabel('Longitudinal velocity u', 'FontSize', 15);
zlabel('\delta(x)', 'FontSize', 20);
view([-146 11])
ylim([0 3])
zlim([0, 2*max(y{j})])


%% normal velocity
figure(10)
v = U_inf./sqrt(2*U_inf*x*rho/mu).*(eta.*y2-y1);
% at x = 0.1
y = eta*sqrt(0.1)/sqrt(rho*U_inf/2/mu);
plot(v, y, 'LineWidth', 2, 'color', 'black')
title(' Normal Velocity of Flat Plate at x = 0.1 ', 'FontSize', 14);
xlabel('V', 'FontSize', 20);
ylabel('Y', 'FontSize', 20);
axis tight
grid on