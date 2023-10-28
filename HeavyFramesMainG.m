clear
clc
close all

R = 1.393; % radius
V = 132174; % vertical load per U/C strut previously 147445
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

% right gear
alpha = acos(y/R); % angle of attachment of gear strut
x = sqrt(R^2 - y^2); % cabin floor half width
Pl = V*sin(alpha); % tangential component of gear load
Ql = V*cos(alpha); % radial component of gear load
Tl = (WS-x)*V; % resultant torque at gear attachment, +ve anticlockwise

Ml1 = Ql*R*W1 + Pl*R*W4 + Tl*W7;
Nl1 = Ql*W2 + Pl*W5 + Tl*W8/R;
Sl1 = Ql*W3 + Pl*W6 + Tl*W9/R;


% left gear
alpha = acos(y/R); % angle of attachment of gear strut
x = sqrt(R^2 - y^2); % cabin floor half width
P2 = -V*sin(alpha); % tangential component of gear load
Q2 = V*cos(alpha); % radial component of gear load
T2 = -(WS-x)*V; % resultant torque at gear attachment, +ve anticlockwise

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


shiftNum = round(shiftNum/2);
MlPlot(1:shiftNum) = MlTot(1001-shiftNum:1000);
MlPlot(shiftNum+1:1000) = MlTot(1:1000-shiftNum);
NlPlot(1:shiftNum) = NlTot(1001-shiftNum:1000);
NlPlot(shiftNum+1:1000) = NlTot(1:1000-shiftNum);
SlPlot(1:shiftNum) = SlTot(1001-shiftNum:1000);
SlPlot(shiftNum+1:1000) = SlTot(1:1000-shiftNum);

figure
plot(phi,MlPlot)
hold on
plot(phi,NlPlot)
plot(phi,SlPlot)
hold off
legend('M','N','S')

figure
polarplot(linspace(0,2*pi,1000),ones(1,1000)*2*10^5)
hold on
polarplot([phi,2*pi]+alpha-pi/2,[MlTot,MlTot(1)]+2*10^5)
hold off
set(gca,'rTick', [])
thetaticklabels({'270^{\circ}','300^{\circ}','330^{\circ}','0^{\circ}','30^{\circ}','60^{\circ}','90^{\circ}','120^{\circ}','150^{\circ}','180^{\circ}','210^{\circ}','240^{\circ}'});
legend('U/C Heavy Frame','Bending Moment','location','northeastoutside')

Mmax = max(abs(MlTot))
Nmax = max(abs(NlTot))
Smax = max(abs(SlTot))





