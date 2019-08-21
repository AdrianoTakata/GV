function [f,y,z] = Runge_Kutta4(a)

f1 = @(x, f, y, z) y;
f2 = @(x, f, y, z) z;
f3 = @(x, f, y, z) -f*z;

f(1) = 0.0;
y(1) = 0.0;
z(1) = a;

h = 0.01;
eta = 0:h:8.8;

for i = 1:(length(eta)-1)
a = h.*[f1(eta(i), f(i), y(i), z(i)), f2(eta(i), f(i), y(i), z(i)), f3(eta(i), f(i), y(i), z(i))];
b = h.*[f1(eta(i)+h/2, f(i)+a(1)/2, y(i)+a(2)/2, z(i)+a(3)/2), f2(eta(i)+h/2, f(i)+a(1)/2, y(i)+a(2)/2, z(i)+a(3)/2), f3(eta(i)+h/2, f(i)+a(1)/2, y(i)+a(2)/2, z(i)+a(3)/2)];
c = h.*[f1(eta(i)+h/2, f(i)+b(1)/2, y(i)+b(2)/2, z(i)+b(3)/2), f2(eta(i)+h/2, f(i)+b(1)/2, y(i)+b(2)/2, z(i)+b(3)/2), f3(eta(i)+h/2, f(i)+b(1)/2, y(i)+b(2)/2, z(i)+b(3)/2)];
d = h.*[f1(eta(i)+h, f(i)+c(1), y(i)+c(2), z(i)+c(3)), f2(eta(i)+h, f(i)+c(1), y(i)+c(2), z(i)+c(3)), f3(eta(i)+h, f(i)+c(1), y(i)+c(2), z(i)+c(3))];
z(i+1) = z(i)+ 1/6*(a(3)+2*b(3)+2*c(3)+d(3));
y(i+1) = y(i)+ 1/6*(a(2)+2*b(2)+2*c(2)+d(2));
f(i+1) = f(i)+ 1/6*(a(1)+2*b(1)+2*c(1)+d(1));
end