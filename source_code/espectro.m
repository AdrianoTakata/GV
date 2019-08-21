clc;
close all;
clear all;
format long;

sonda01=load('sonda_l3_02.dat');

ini = 20000;
[fin,c] = size(sonda01);
fin = 40101;

sonda01 = sonda01(ini:fin,:);
ini = 1;
[fin,c] = size(sonda01);

x1 = hanning(fin);
size(x1)

x2 = x1.*sonda01(:,2);
size(x2)

%%x3 = 2.*(fft(x2))./(fin/2);
x3 = fft(x2);
x3 = x3.*conj(x3)/fin;
size(x3)
x3 = x3(1:fin/2);
x3(1:2) = [0; 0];

t  = sonda01(fin,1) - sonda01(1,1)
dt = sonda01(2,1) - sonda01(1,1)

%%f = 0:1/t:1/0.0001;
faq = 1/dt;
f = (0:fin/2) * faq / fin;
size(f)

loglog(f(1:fin/2), x3(1:fin/2), '-', f(1:fin/2), f(1:fin/2).^(-5/3), '-');
ylabel('E(f)'), xlabel('Frequencia'), title('Espectro de Energia');

f_pot = f.^(-5/3)';
M = [f(1:end-1)', x3, f_pot(1:end-1)];
M = M(2:end, :);
save('espectro_l3_02.dat', 'M', '-ascii', '-double');





