% HTP frame

clear
clc
clear all

R = 1.266/2; % radius
Q = 43700/2;
P = 23650;
T = -28909;


phi = linspace(0,2*pi,1000);

W1 = ((cos(phi)/2)-(pi-phi).*sin(phi)+1)/(2*pi);
W2 = ((3*cos(phi)/2)+(pi-phi).*sin(phi))/(2*pi);
W3 = ((pi-phi).*cos(phi)-(sin(phi)/2))/(2*pi);

W4 = (1 / (2 * pi)) * (3*sin(phi)/2 + (pi - phi).*(cos(phi) - 1));
W5 = (1 / (2 * pi)) * (sin(phi)/2 - (pi-phi).*cos(phi));
W6 = (1 / (2 * pi)) * ((pi - phi).* sin(phi) - 1 - 0.5 * cos(phi));

W7 = (pi-2*sin(phi)-phi)/(2*pi);
%W8 = ((3*cos(phi)/2)+(pi-phi).*sin(phi))/(2*pi);
W8 = sin(phi)/pi;
W9 = (1 / (2 * pi)) * (1 + 2 * cos(phi));

M = Q*R*W1 + P*R*W4 + T*W7;
N = Q*W2 + P*W5 + T*W8/R;
S = Q*W3 + P*W6 + T*W9/R;


figure
plot(phi,M)
hold on
plot(phi,N)
plot(phi,S)
hold off
legend('M','N','S')

Mmax = max(abs(M))
Nmax = max(abs(N))
Smax = max(abs(S))