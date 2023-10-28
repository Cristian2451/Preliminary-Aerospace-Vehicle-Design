% Code to optimise a fuselage structure for a given material,
% skin thickness, stringer separation and max stringer direct stress

clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters set in conceptual design:
R = 1.393; % fuselage outer radius
Fa = -732665; % max axial load on fuselage (wing spar)
Ft = 0; % max tangential load on fuselage
T = 26424; % max torque experienced by fuselage
BM = -2.50E+06; % max bending moment in fuselage

% material properties: Al-2090-T83
%SigY = 2.31e8;
SigY = 4.84e8; % tensile yield stress
SigShearY = SigY/sqrt(3); % shear yield stress
%E = 7.385e10;
E = 7.79e10; % Young's modulus

% parameters:
n = 400; % number of values to test for each parameter
tMin = 1e-3; % minimum skin thickness to test
tMax = 3e-3;
SigMin = SigY/5; % minimum stringer direct stress to test
SigMax = SigY/2;
NsMin = floor(2*pi*R/0.5); % minimum number of stringers to test
NsMax = ceil(2*pi*R/0.08); % 0.08

% frame parameters:
Tframe = 2.5e-3; % frame thickness
Hframe = 5e-2; % frame height
Cf = 1/16000; % empirical constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check minimum skin thickness required for shear stress
phi = linspace(0,360,37);
shear = Fa*sind(phi)/(pi*R);
pointLoad = Ft*cosd(phi)/(pi*R);
torqueLoad = (T+Ft*R)/(2*pi*R^2);
shearTot = shear + pointLoad + torqueLoad;
% find value with max magnitude
shearMax = max([max(shearTot),abs(min(shearTot))]);

% adjust minimum skin thickness to test if necessary
if tMin < shearMax/SigShearY
    tMin = ceil(10000*shearMax/SigShearY)/10000;
end


% determine values of three main parameters to test
sigma = linspace(SigMin,SigMax,n); % max stringer direct stress
Ns = NsMin:1:NsMax; % stringer separation
tSkin = tMin:0.0001:tMax; % skin thickness

l = 1;
m = 1;
result(n^3,7) = 0;
Ds(length(Ns)) = 0;
% iterate over parameters
for i = 1:n % sigma
    disp(['i = ',num2str(i)]) % indicates progress
    for j = 1:length(tSkin) % tSkin
        for k = 1:length(Ns) % Ns
            Ds(k) = 2*pi*R/Ns(k);
            sigCrit = 3.62*E*(tSkin(j)/Ds(k))^2; % critical panel buckling stress
            if sigCrit > sigma(i) % check critical buckling stress excedes sigma
                % stringers
                thetaString = linspace(0,2*pi,Ns(k)+1); % angular distribution of stringers
                thetaString(end)=[]; % remove duplication of zero (2*pi)
                YstringSq = (R*sin(thetaString)).^2; % stringer y positions squared
                FusIreq = R*abs(BM)/sigma(i); % required 2nd moment of area of fuselage
                Astring = FusIreq/sum(YstringSq); % required area per stringer
                Hstring = sqrt(Astring); % stringer height
                Istring = (Hstring^4)/12; % 2nd moment of area of stringer
                Df = sqrt((pi^2)*E*Istring/(Astring*sigma(i))); % frame spacing
                
                if Df > 0.35
                    % light frames
                    Iframe = (Cf*abs(BM)*(2*R)^2)/(Df*E); % required 2nd moment of area of frame
                    Bframe = (Iframe-(Tframe*(Hframe-2*Tframe)^3)/12)/(2*((Tframe^3)/12+Tframe*(Hframe/2-Tframe/2)^2)); % required frame base
                    Aframe = (2*Bframe*Tframe)+Tframe*(Hframe-2*Tframe);
                    
                    % check overall volume and add parameters to array
                    result(l,1) = sigma(i); % max stress in stringers
                    result(l,2) = tSkin(j); % skin thickness
                    result(l,3) = Ns(k); % number of stringers
                    result(l,4) = Hstring; % stringer height & width
                    result(l,5) = Df; % frame spacing
                    result(l,6) = Bframe; % frame base
                    result(l,7) = 2*pi*R*tSkin(j) + Ns(k)*Astring + 2*pi*R*Aframe/Df; % volume per length
                    l = l + 1;
    
                    % values for 3D plot with constant stringer separation
                    if Ns(k) == NsMax
                    %if true
                        plotX(i) = sigma(i);
                        plotY(j) = tSkin(j);
                        plotZ(j,i) = result(l-1,7);
                        %D(m) = result(l-1,7);
                        m = m + 1;
                    end
                end
            end
        end
    end
end

% find minimum mass combination
[minVol,minVolIndex] = min(result(1:l-1,7));
disp(['sigma = ',num2str(result(minVolIndex,1))])
disp(['tSkin = ',num2str(result(minVolIndex,2))])
disp(['Ns = ',num2str(result(minVolIndex,3))])
disp(['Stringer height & width = ',num2str(result(minVolIndex,4))])
disp(['Frame spacing = ',num2str(result(minVolIndex,5))])
disp(['Frame base = ',num2str(result(minVolIndex,6))])
disp(['vol/L = ',num2str(result(minVolIndex,7))])

% plot mass for each iteration
figure
plot(result(1:l-1,7))

% 3D plot with constant stringer separation
figure
plotZ(plotZ==0) = nan;
surf(plotX/(10^6),plotY*1000,plotZ,'EdgeColor','none')  % 'EdgeColor','none','FaceAlpha',0.5
colormap('jet');
%scatter3(plotX/(10^6),plotY*1000,plotZ,'LineWidth',D)
xlabel('Maximum Stringer Direct Stress / MPa')
ylabel('Skin Thickness / mm')
zlabel('Volume per Length / m^2')

