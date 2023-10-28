% plots h vs A for U/C heavy frame

clear
clc
close all

n = 1000;

I = 1.45e-5;
t = 0.02;
h = linspace(0.05,0.1,n);

A(n) = 0;

for i = 1:n
    b = (I - (t*(h(i)-2*t)^3)/12) / (2*(t^3/12)+t*(h(i)/2 - t/2)^2);
    A(i) = 2*b*t + (h(i)-2*t)*t;
end

plot(h,A)
xlabel('Frame Cross-Section Height / m')
ylabel('Frame Cross-Sectional Area / m^2')