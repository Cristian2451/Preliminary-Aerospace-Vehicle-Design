clear
clc


R = 1.393; % radius
V = 147445; % vertical load per U/C strut
WS = 2.11; % distance from centre line to main gear wheels
y = 0.653; % distance from centre of a/c to cabin floor


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

alpha = acos(y/R); % angle of attachment of gear strut
x = sqrt(R^2 - y^2); % cabin floor half width
Pl = V*sin(alpha); % tangential component of gear load
Ql = V*cos(alpha); % radial component of gear load
Tl = (WS-x)*V; % resultant torque at gear attachment, +ve anticlockwise

Ml1 = Ql*R*W1 + Pl*R*W4 + Tl*W7;
Nl1 = Ql*W2 + Pl*W5 + Tl*W8/R;
Sl1 = Ql*W3 + Pl*W6 + Tl*W9/R;

% redefine phi to account for difference in position of two struts
phi = linspace(2*alpha,2*pi+2*alpha,1000);

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

Ml2 = Ql*R*W1 - Pl*R*W4 - Tl*W7;
Nl2 = Ql*W2 - Pl*W5 - Tl*W8/R;
Sl2 = Ql*W3 - Pl*W6 - Tl*W9/R;

MlTot = Ml1 + Ml2;
NlTot = Nl1 + Nl2;
SlTot = Sl1 + Sl2;

phi = linspace(0,2*pi,1000);


figure
plot(phi,MlTot)
hold on
plot(phi,NlTot)
plot(phi,SlTot)
hold off
legend('M','N','S')

Mmax = max(abs(MlTot))
Nmax = max(abs(NlTot))
Smax = max(abs(SlTot))





