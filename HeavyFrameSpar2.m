% rear wing spar

clear
clc
clear all

R = 1.393; % radius
V = -732665/2; % vertical load per wing strut
alpha = 20*pi/180;


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

% right attachment
Pl = V*sin(alpha); % tangential component of wing load
Ql = V*cos(alpha); % radial component of wing load
Tl = 0; % resultant torque at attachment, +ve anticlockwise

Ml1 = Ql*R*W1 + Pl*R*W4 + Tl*W7;
Nl1 = Ql*W2 + Pl*W5 + Tl*W8/R;
Sl1 = Ql*W3 + Pl*W6 + Tl*W9/R;


% left attachment
P2 = -V*sin(alpha); % tangential component of gear load
Q2 = V*cos(alpha); % radial component of gear load
T2 = 0; % resultant torque at attachment, +ve anticlockwise

Ml2 = Q2*R*W1 + P2*R*W4 + T2*W7;
Nl2 = Q2*W2 + P2*W5 + T2*W8/R;
Sl2 = Q2*W3 + P2*W6 + T2*W9/R;


shiftNum = round(1000*(2*alpha/(2*pi)));

MlTot(1:1000-shiftNum) = Ml1(1:1000-shiftNum) + Ml2(shiftNum+1:1000);
MlTot(1001-shiftNum:1000) = Ml1(1001-shiftNum:1000) + Ml2(1:shiftNum);

NlTot(1:1000-shiftNum) = Nl1(1:1000-shiftNum) + Nl2(shiftNum+1:1000);
NlTot(1001-shiftNum:1000) = Nl1(1001-shiftNum:1000) + Nl2(1:shiftNum);

SlTot(1:1000-shiftNum) = Sl1(1:1000-shiftNum) + Sl2(shiftNum+1:1000);
SlTot(1001-shiftNum:1000) = Sl1(1001-shiftNum:1000) + Sl2(1:shiftNum);


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