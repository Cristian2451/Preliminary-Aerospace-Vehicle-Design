clear
clc

x_FA = 0.375; % Flexural Axis 
x_AC = 0.25; % Aerodynamic Centre
x_CG = 0.4; % Wing + Fuel Centre of Gravity (update to more accurate one!!)

MTOW = 385434.9; % [N]
W_NFW = MTOW - 7475*9.81; % No fuel at wing
Sref = 79.471;
b = 28.2;
W_S = 4850;
c_bar = 3.25;
sweep_c4 = 12;
taper = 0.4;
rho = 1.225;
AR = 10;
c_root = 2*Sref/(b*(1+taper));
c_tip = c_root*taper;
sweep_FA = atand(tand(sweep_c4)-(4/AR)*((x_FA - x_AC)*(1-taper)/(1+taper)));
sweep_LE = atand(tand(sweep_c4)-(4/AR)*((0 - x_AC)*(1-taper)/(1+taper)));
sweep_hinge = atand(tand(sweep_c4)-(4/AR)*((0.75 - x_AC)*(1-taper)/(1+taper)));
sweep_10 = atand(tand(sweep_c4)-(4/AR)*((0.1 - x_AC)*(1-taper)/(1+taper)));

n = 3.75;
n_A = 2.5;
n_A_low = -1;
n_low = -1.5;
s = b*0.5; %/cosd(sweep_c4); % Should this be at c/4 or flexural axis?????
y = linspace(0,s,100);


Clmax_clean = 1.365;
Clmin_airfoil = -0.8;
Clmin = 0.6*Clmin_airfoil*cosd(sweep_c4);
Cl_alpha = 5.84;


V_cruise = 143; % [m/s] for first leg (greates EAS)
V_ne = (1/0.8)*(V_cruise/0.69)*(0.69+0.07); %formula from CS_25 





% x_airfoil = [1, 0.95033, 0.90066, 0.8509, 0.80103, 0.75107, 0.70101, 0.65086, 0.60064, 0.55035, 0.5, 0.44962, 0.39923, 0.34884, 0.29846, 0.24881, 0.19781, 0.14757, 0.09746, 0.07247, 0.04757, 0.02283, 0.01059, 0.0058, 0.00347, 0, 0.00653, 0.0092, 0.01441, 0.02717, 0.05243, 0.07753, 0.10254, 0.15243, 0.20219, 0.25189, 0.30154, 0.35116, 0.40077, 0.45038, 0.5, 0.54965, 0.59936, 0.64914, 0.69899, 0.74893, 0.79897, 0.8491, 0.89934, 0.94967, 1];
% y_airfoil = [0, 0.00986, 0.01979, 0.02974, 0.03935, 0.04847, 0.05686, 0.0644, 0.07085, 0.07602, 0.07963, 0.08139, 0.08139, 0.07971, 0.07658, 0.07193, 0.06562, 0.05741, 0.04672, 0.0401, 0.03227, 0.02234, 0.01588, 0.01236, 0.0101, 0, -0.0081, -0.00956, -0.0116, -0.0149, -0.01963, -0.02314, -0.02604, -0.03049, -0.03378, -0.03613, -0.0377, -0.03851, -0.03855, -0.03759, -0.03551, -0.03222, -0.02801, -0.0232, -0.01798, -0.01267, -0.00751, -0.00282, 0.00089, 0.00278, 0];
% 
% figure 
% plot([x_FS,x_FS,x_RS,x_RS,x_FS],[y_avg+0.02,-y_avg+0.02,-y_avg+0.02,y_avg+0.02,y_avg+0.02],'b-','LineWidth',2)
% hold on
% plot(x_FA,0.02,'b*','LineWidth',2,'MarkerSize',14)
% %plot(x_CG,0,'ko','LineWidth',2,'MarkerSize',14)
% plot(x_airfoil,y_airfoil,'k-','LineWidth',2)
% hold off
% grid minor
% box on
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',18)
% set(gca,'TickLabelInterpreter','latex')
% xlabel('x/c [m]','Interpreter', 'latex', 'FontSize',24)
% ylabel('y/c [m]','Interpreter', 'latex', 'FontSize',24)
% legend('Wingbox','Flexural Axis','Interpreter', 'latex', 'FontSize',22)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Limit Load Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_up = 2.5;
n_down = -1;
V_A_up = sqrt(MTOW*n_up/(0.5*rho*Clmax_clean*Sref));
V_A_down = sqrt(MTOW*n_down/(0.5*rho*Clmin*Sref));

V_left_up = linspace(0,V_A_up,10);
V_left_down = linspace(0,V_A_down,10);
n_left_up = 0.5*rho*V_left_up.^2*Clmax_clean*Sref/MTOW;
n_left_down = 0.5*rho*V_left_down.^2*Clmin*Sref/MTOW;


%%%%%%%%%%%%%%%%%%%%%%%%%% Ultimate Load Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_up_UL = 2.5*1.5;
n_down_UL = -1*1.5;
V_A_up_UL = sqrt(MTOW*n_up_UL/(0.5*rho*Clmax_clean*Sref));
V_A_down_UL = sqrt(MTOW*n_down_UL/(0.5*rho*Clmin*Sref));

V_left_up_UL = linspace(0,V_A_up_UL,30);
V_left_down_UL = linspace(0,V_A_down_UL,30);
n_left_up_UL = 0.5*rho*V_left_up_UL.^2*Clmax_clean*Sref/MTOW;
n_left_down_UL = 0.5*rho*V_left_down_UL.^2*Clmin*Sref/MTOW;


figure
plot(V_left_up_UL,n_left_up_UL,'b-','LineWidth',1.5)
hold on
plot(V_left_up,n_left_up,'k--','LineWidth',1.5)
plot(V_left_down_UL,n_left_down_UL,'b-','LineWidth',1.5)
plot([V_A_up_UL,V_ne,V_ne,V_A_down_UL],[n_up_UL,n_up_UL,n_down_UL,n_down_UL],'b-','LineWidth',1.5)
plot(V_left_down,n_left_down,'k--','LineWidth',1.5)
plot([V_A_up,V_ne,V_ne,V_cruise,V_A_down],[n_up,n_up,0,n_down,n_down],'k--','LineWidth',1.5)
p1=plot(120.426,2.5,'ko',MarkerFaceColor='k');
label(p1,'$V_A$','Interpreter', 'latex','FontSize',15)
p2=plot(196.88,3.75,'ko',MarkerFaceColor='k');
label(p2,'$V_D$','Interpreter', 'latex','FontSize',15)
hold off
grid on 
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('EAS [m/s]','Interpreter', 'latex', 'FontSize',18)
ylabel('Load Factor n','Interpreter', 'latex', 'FontSize',18)
legend('Ultimate Load','Limit Load','Interpreter', 'latex', 'FontSize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gust Loads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = 1.225;
mu = W_S/(0.5*rho*9.81*c_bar*Cl_alpha);
K = 0.88*mu/(5.3 + mu);

U_de_turb = 20;
U_de_cruise = 15.2;
U_de_dive = 7.6;

U_turb = K*U_de_turb;
U_cruise = K*U_de_cruise;
U_dive = K*U_de_dive;

V_s1 = sqrt(MTOW*1/(0.5*rho*Clmax_clean*Sref)); % Stall speed at level flight

V_g = V_s1*sqrt(1 + (K*U_de_turb*V_cruise*1.94384*Cl_alpha)/(498*W_S*0.204816));

% dn_turb = 0.5*rho*U_turb*V_g*Cl_alpha/W_S;
% dn_cruise = 0.5*rho*U_cruise*V_cruise*Cl_alpha/W_S;
% dn_dive = 0.5*rho*U_dive*V_ne*Cl_alpha/W_S;

dn_turb = Cl_alpha*0.5*rho*V_g*U_turb/W_S;
dn_cruise = 0.5*rho*U_cruise*V_cruise*Cl_alpha/W_S;
dn_dive = 0.5*rho*U_dive*V_ne*Cl_alpha/W_S;

figure
plot(V_left_up_UL,n_left_up_UL,'b-','LineWidth',1.5)
hold on
plot(V_left_up,n_left_up,'k--','LineWidth',1.5)
plot([0,V_g,V_cruise,V_ne,V_ne,V_cruise,V_g,0],[1,1+dn_turb,1+dn_cruise,1+dn_dive,1-dn_dive,1-dn_cruise,1-dn_turb,1],'r-.','LineWidth',1.5)
plot(V_left_down_UL,n_left_down_UL,'b-','LineWidth',1.5)
plot([V_A_up_UL,V_ne,V_ne,V_A_down_UL],[n_up_UL,n_up_UL,n_down_UL,n_down_UL],'b-','LineWidth',1.5)
plot(V_left_down,n_left_down,'k--','LineWidth',1.5)
plot([V_A_up,V_ne,V_ne,V_cruise,V_A_down],[n_up,n_up,0,n_down,n_down],'k--','LineWidth',1.5)
p1=plot(120.426,2.5,'ko',MarkerFaceColor='k');
label(p1,'$V_A$','Interpreter', 'latex','FontSize',15)
p2=plot(196.88,3.75,'ko',MarkerFaceColor='k');
label(p2,'$V_D$','Interpreter', 'latex','FontSize',15)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('EAS [m/s]','Interpreter', 'latex', 'FontSize',18)
ylabel('Load Factor n','Interpreter', 'latex', 'FontSize',18)
legend('Ultimate Load','Limit Load','Gust Load','Interpreter', 'latex', 'FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%% Normal Landing Load %%%%%%%%%%%%%%%%%%%%%%%%%%%

V_v = 3.05;
n_gear = 3;
W_land = 356806.6;
KE_v = 0.5*(W_land/9.81)*V_v^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inertia Loads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3.5 Load Factor

W_fuel = zeros(1,length(y));
W_wing = zeros(1,length(y));
W_fuel_dist = zeros(1,length(y));
W_wing_dist = zeros(1,length(y));
y_n = linspace(0,b/2,100);
W_engine = zeros(1,length(y));
W_engine(length(y)*0.3) = 14601*n;
W_dist_engine = zeros(1,length(y));
W_dist_engine(length(y)*0.3) = 14601;
Thrust_dist = zeros(1,length(y));
Thrust_dist(length(y)*0.3) = 63000;

for i = 1:length(y_n)
    if y_n(i) < 13
        W_fuel(i) = (3737.5 - (-22.368*y_n(i)^2 + 572.56*y_n(i) - 6.2382))*9.81*n;
    end
    W_wing(i) = (2184 - (-8.7978*y_n(i)^2 + 277.46*y_n(i) + 21.242))*9.81*n;
end

for i = 1:length(y_n)
    if y_n(i) < 13
        W_fuel_dist(i) = (3737.5 - (-22.368*y_n(i)^2 + 572.56*y_n(i) - 6.2382))*9.81;
    end
    W_wing_dist(i) = (2184 - (-8.7978*y_n(i)^2 + 277.46*y_n(i) + 21.242))*9.81;
end

W_fuel_dist(1:92) = W_fuel_dist(1:92) - 788;
for i = 1:length(y_n)-1
    W_wing_dist(length(y_n) - i) = W_wing_dist(length(y_n) - i) - W_wing_dist(length(y_n) - (i-1));
    W_fuel_dist(length(y_n) - i) = W_fuel_dist(length(y_n) - i) - W_fuel_dist(length(y_n) - (i-1));
end

SF = zeros(1,length(y));
dM = zeros(1,length(y));
BM = zeros(1,length(y));

SF_ZF = zeros(1,length(y));
dM_ZF = zeros(1,length(y));
BM_ZF = zeros(1,length(y));

for i = 1:length(y)-1
    SF(i) = W_wing(i) + W_fuel(i) + sum(W_engine(i:end));
    SF_ZF(i) = W_wing(i) + sum(W_engine(i:end));
end

for i = 1:length(y)-1
    dM(i) = (SF(i+1) + SF(i))*(y(i+1) - y(i))*0.5;
    dM_ZF(i) = (SF_ZF(i+1) + SF_ZF(i))*(y(i+1) - y(i))*0.5;
end

for i = 1:length(y)-1
    BM(i) = sum(dM(i:end));
    BM_ZF(i) = sum(dM_ZF(i:end));
end


% 2.5 Load Factor (V_A)

W_fuel_A = zeros(1,length(y));
W_wing_A = zeros(1,length(y));
y_n = linspace(0,b/2,100);
W_engine_A = zeros(1,length(y));
W_engine_A(length(y)*0.3) = 14601*n_A;

for i = 1:length(y_n)
    if y_n(i) < 13
        W_fuel_A(i) = (3737.5 - (-22.368*y_n(i)^2 + 572.56*y_n(i) - 6.2382))*9.81*n_A;
    end
    W_wing_A(i) = (2184 - (-8.7978*y_n(i)^2 + 277.46*y_n(i) + 21.242))*9.81*n_A;
end

SF_A = zeros(1,length(y));
dM_A = zeros(1,length(y));
BM_A = zeros(1,length(y));

SF_ZF_A = zeros(1,length(y));
dM_ZF_A = zeros(1,length(y));
BM_ZF_A = zeros(1,length(y));

for i = 1:length(y)-1
    SF_A(i) = W_wing_A(i) + W_fuel_A(i) + sum(W_engine_A(i:end));
    SF_ZF_A(i) = W_wing_A(i) + sum(W_engine_A(i:end));
end

for i = 1:length(y)-1
    dM_A(i) = (SF_A(i+1) + SF_A(i))*(y(i+1) - y(i))*0.5;
    dM_ZF_A(i) = (SF_ZF_A(i+1) + SF_ZF_A(i))*(y(i+1) - y(i))*0.5;
end

for i = 1:length(y)-1
    BM_A(i) = sum(dM_A(i:end));
    BM_ZF_A(i) = sum(dM_ZF_A(i:end));
end

% -1.5 Load Factor

W_fuel_low = zeros(1,length(y));
W_wing_low = zeros(1,length(y));
y_n = linspace(0,b/2,100);
W_engine_low = zeros(1,length(y));
W_engine_low(length(y)*0.3) = 14601*n_low;

for i = 1:length(y_n)
    if y_n(i) < 13
        W_fuel_low(i) = (3737.5 - (-22.368*y_n(i)^2 + 572.56*y_n(i) - 6.2382))*9.81*n_low;
    end
    W_wing_low(i) = (2184 - (-8.7978*y_n(i)^2 + 277.46*y_n(i) + 21.242))*9.81*n_low;
end

SF_low = zeros(1,length(y));
dM_low = zeros(1,length(y));
BM_low = zeros(1,length(y));

SF_ZF_low = zeros(1,length(y));
dM_ZF_low = zeros(1,length(y));
BM_ZF_low = zeros(1,length(y));

for i = 1:length(y)-1
    SF_low(i) = W_wing_low(i) + W_fuel_low(i) + sum(W_engine_low(i:end));
    SF_ZF_low(i) = W_wing_low(i) + sum(W_engine_low(i:end));
end

for i = 1:length(y)-1
    dM_low(i) = (SF_low(i+1) + SF_low(i))*(y(i+1) - y(i))*0.5;
    dM_ZF_low(i) = (SF_ZF_low(i+1) + SF_ZF_low(i))*(y(i+1) - y(i))*0.5;
end

for i = 1:length(y)-1
    BM_low(i) = sum(dM_low(i:end));
    BM_ZF_low(i) = sum(dM_ZF_low(i:end));
end


% -1 Load Factor (V_A)

W_fuel_A_low = zeros(1,length(y));
W_wing_A_low = zeros(1,length(y));
y_n = linspace(0,b/2,100);
W_engine_A_low = zeros(1,length(y));
W_engine_A_low(length(y)*0.3) = 14601*n_A_low;

for i = 1:length(y_n)
    if y_n(i) < 13
        W_fuel_A_low(i) = (3737.5 - (-22.368*y_n(i)^2 + 572.56*y_n(i) - 6.2382))*9.81*n_A_low;
    end
    W_wing_A_low(i) = (2184 - (-8.7978*y_n(i)^2 + 277.46*y_n(i) + 21.242))*9.81*n_A_low;
end

SF_A_low = zeros(1,length(y));
dM_A_low = zeros(1,length(y));
BM_A_low = zeros(1,length(y));

SF_ZF_A_low = zeros(1,length(y));
dM_ZF_A_low = zeros(1,length(y));
BM_ZF_A_low = zeros(1,length(y));

for i = 1:length(y)-1
    SF_A_low(i) = W_wing_A_low(i) + W_fuel_A_low(i) + sum(W_engine_A_low(i:end));
    SF_ZF_A_low(i) = W_wing_A_low(i) + sum(W_engine_A_low(i:end));
end

for i = 1:length(y)-1
    dM_A_low(i) = (SF_A_low(i+1) + SF_A_low(i))*(y(i+1) - y(i))*0.5;
    dM_ZF_A_low(i) = (SF_ZF_A_low(i+1) + SF_ZF_A_low(i))*(y(i+1) - y(i))*0.5;
end

for i = 1:length(y)-1
    BM_A_low(i) = sum(dM_A_low(i:end));
    BM_ZF_A_low(i) = sum(dM_ZF_A_low(i:end));
end


% figure
% plot(y,SF,'b-','LineWidth',1.5)
% hold on
% plot(y,SF_ZF,'b-.','LineWidth',1.5)
% hold off
% grid on
% box on
% xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',14)
% ylabel('Inertia Shear Force [N]','Interpreter', 'latex', 'FontSize',14)
% legend('Load Case 1 (MTOW)','Load Case 1 (Zero Wing Fuel)','Interpreter', 'latex', 'FontSize',12)
% 
% figure
% plot(y,BM,'b-','LineWidth',1.5)
% hold on
% plot(y,BM_ZF,'b-.','LineWidth',1.5)
% hold off
% grid on
% box on
% xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',14)
% ylabel('Inertia Bending Moment [Nm]','Interpreter', 'latex', 'FontSize',14)
% legend('Load Case 1 (MTOW)','Load Case 1 (Zero Wing Fuel)','Interpreter', 'latex', 'FontSize',12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Aero Loads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_dist = MTOW;
L0_dist = 2*L_dist/(pi*s);
dL_dist = L0_dist*sqrt(1-(y/s).^2);

L = n*MTOW;
L0 = 2*L/(pi*s);
dL = L0*sqrt(1-(y/s).^2);

L_ZF = n*W_NFW;
L0_ZF = 2*L_ZF/(pi*s);
dL_ZF = L0_ZF*sqrt(1-(y/s).^2);

L_A = n_A*MTOW;
L0_A = 2*L_A/(pi*s);
dL_A = L0_A*sqrt(1-(y/s).^2);

L_A_low = n_A_low*MTOW;
L0_A_low = 2*L_A_low/(pi*s);
dL_A_low = L0_A_low*sqrt(1-(y/s).^2);

L_ZF_A = n_A*W_NFW;
L0_ZF_A = 2*L_ZF_A/(pi*s);
dL_ZF_A = L0_ZF_A*sqrt(1-(y/s).^2);

L_low = n_low*MTOW;
L0_low = 2*L_low/(pi*s);
dL_low = L0_low*sqrt(1-(y/s).^2);

L_ZF_low = n_low*W_NFW;
L0_ZF_low = 2*L_ZF_low/(pi*s);
dL_ZF_low = L0_ZF_low*sqrt(1-(y/s).^2);

Lift = zeros(1,length(y));
SF_aero = zeros(1,length(y));
dM_aero = zeros(1,length(y));
BM_aero = zeros(1,length(y));

Lift_ZF = zeros(1,length(y));
SF_aero_ZF = zeros(1,length(y));
dM_aero_ZF = zeros(1,length(y));
BM_aero_ZF = zeros(1,length(y));

Lift_A = zeros(1,length(y));
SF_aero_A = zeros(1,length(y));
dM_aero_A = zeros(1,length(y));
BM_aero_A = zeros(1,length(y));

Lift_A_low = zeros(1,length(y));
SF_aero_A_low = zeros(1,length(y));
dM_aero_A_low = zeros(1,length(y));
BM_aero_A_low = zeros(1,length(y));

Lift_ZF_A = zeros(1,length(y));
SF_aero_ZF_A = zeros(1,length(y));
dM_aero_ZF_A = zeros(1,length(y));
BM_aero_ZF_A = zeros(1,length(y));

Lift_low = zeros(1,length(y));
SF_aero_low = zeros(1,length(y));
dM_aero_low = zeros(1,length(y));
BM_aero_low = zeros(1,length(y));

Lift_ZF_low = zeros(1,length(y));
SF_aero_ZF_low = zeros(1,length(y));
dM_aero_ZF_low = zeros(1,length(y));
BM_aero_ZF_low = zeros(1,length(y));

for i = 1:length(y)-1
    Lift(i) = (dL(i+1) + dL(i))*(y(i+1) - y(i))*0.5;
    Lift_ZF(i) = (dL_ZF(i+1) + dL_ZF(i))*(y(i+1) - y(i))*0.5;
    Lift_A(i) = (dL_A(i+1) + dL_A(i))*(y(i+1) - y(i))*0.5;
    Lift_A_low(i) = (dL_A_low(i+1) + dL_A_low(i))*(y(i+1) - y(i))*0.5;
    Lift_ZF_A(i) = (dL_ZF_A(i+1) + dL_ZF_A(i))*(y(i+1) - y(i))*0.5;
    Lift_low(i) = (dL_low(i+1) + dL_low(i))*(y(i+1) - y(i))*0.5;
    Lift_ZF_low(i) = (dL_ZF_low(i+1) + dL_ZF_low(i))*(y(i+1) - y(i))*0.5;
end

for i = 1:length(y)-1
    SF_aero(i) = sum(Lift(i:end));
    SF_aero_ZF(i) = sum(Lift_ZF(i:end));
    SF_aero_A(i) = sum(Lift_A(i:end));
    SF_aero_A_low(i) = sum(Lift_A_low(i:end));
    SF_aero_ZF_A(i) = sum(Lift_ZF_A(i:end));
    SF_aero_low(i) = sum(Lift_low(i:end));
    SF_aero_ZF_low(i) = sum(Lift_ZF_low(i:end));
end

for i = 1:length(y)-1
    dM_aero(i) = (SF_aero(i+1) + SF_aero(i))*(y(i+1) - y(i))*0.5;
    dM_aero_ZF(i) = (SF_aero_ZF(i+1) + SF_aero_ZF(i))*(y(i+1) - y(i))*0.5;
    dM_aero_A(i) = (SF_aero_A(i+1) + SF_aero_A(i))*(y(i+1) - y(i))*0.5;
    dM_aero_A_low(i) = (SF_aero_A_low(i+1) + SF_aero_A_low(i))*(y(i+1) - y(i))*0.5;
    dM_aero_ZF_A(i) = (SF_aero_ZF_A(i+1) + SF_aero_ZF_A(i))*(y(i+1) - y(i))*0.5;
    dM_aero_low(i) = (SF_aero_low(i+1) + SF_aero_low(i))*(y(i+1) - y(i))*0.5;
    dM_aero_ZF_low(i) = (SF_aero_ZF_low(i+1) + SF_aero_ZF_low(i))*(y(i+1) - y(i))*0.5;
end

for i = 1:length(y)-1
    BM_aero(i) = sum(dM_aero(i:end));
    BM_aero_ZF(i) = sum(dM_aero_ZF(i:end));
    BM_aero_A(i) = sum(dM_aero_A(i:end));
    BM_aero_A_low(i) = sum(dM_aero_A_low(i:end));
    BM_aero_ZF_A(i) = sum(dM_aero_ZF_A(i:end));
    BM_aero_low(i) = sum(dM_aero_low(i:end));
    BM_aero_ZF_low(i) = sum(dM_aero_ZF_low(i:end));
end

% figure
% plot(y,SF_aero,'b-','LineWidth',1.5)
% hold on
% plot(y,SF_aero_ZF,'b-.','LineWidth',1.5)
% hold off
% grid on
% box on
% xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',14)
% ylabel('Aerodynamic Shear Force [N]','Interpreter', 'latex', 'FontSize',14)
% legend('Load Case 1 (MTOW)','Load Case 1 (Zero Wing Fuel)','Interpreter', 'latex', 'FontSize',12)
% 
% figure
% plot(y,BM_aero,'b-','LineWidth',1.5)
% hold on
% plot(y,BM_aero_ZF,'b-.','LineWidth',1.5)
% hold off
% grid on
% box on
% xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',14)
% ylabel('Aerodynamic Bending Moment [Nm]','Interpreter', 'latex', 'FontSize',14)
% legend('Load Case 1 (MTOW)','Load Case 1 (Zero Wing Fuel)','Interpreter', 'latex', 'FontSize',12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Torque  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Moment defined as +ve in anti-clockwise direction

Cm0 = -0.023;
c = (c_tip - c_root)*(y/s) + c_root;

dM0 = 0.5*rho*V_ne^2*Sref*2.991*Cm0*Lift/(s*n*MTOW/2); % Moment per unit span (correction for sweep?)
dM0_A = 0.5*rho*V_A_up^2*Sref*2.991*Cm0*Lift/(s*n_A*MTOW/2); % Moment per unit span (correction for sweep?)
dM0_A_low = 0.5*rho*V_A_down^2*Sref*2.991*Cm0*Lift/(s*n_A_low*MTOW/2); % Moment per unit span (correction for sweep?)
dM0_low = 0.5*rho*V_ne^2*Sref*2.991*Cm0*Lift/(s*n*MTOW/2); % Moment per unit span (correction for sweep?)

dM_lift = -dL.*c*(x_FA - x_AC);
dM_lift_low = -dL_low.*c*(x_FA - x_AC);
dM_lift_A = -dL_A.*c*(x_FA - x_AC);
dM_lift_A_low = -dL_A_low.*c*(x_FA - x_AC);
dM_lift_ZF = -dL_ZF.*c*(x_FA - x_AC);
dM_lift_ZF_low = -dL_ZF_low.*c*(x_FA - x_AC);

M0 = zeros(1,length(y));
M0_A = zeros(1,length(y));
M0_A_low = zeros(1,length(y));
M0_low = zeros(1,length(y));
M_lift = zeros(1,length(y));
M_lift_low = zeros(1,length(y));
M_lift_A = zeros(1,length(y));
M_lift_A_low = zeros(1,length(y));
M_lift_ZF = zeros(1,length(y));
M_lift_ZF_low = zeros(1,length(y));

for i = 1:length(y)-1
    M0(i) = (dM0(i+1) + dM0(i))*(y(i+1) - y(i))*0.5;
    M0_A(i) = (dM0_A(i+1) + dM0_A(i))*(y(i+1) - y(i))*0.5;
    M0_A_low(i) = (dM0_A_low(i+1) + dM0_A_low(i))*(y(i+1) - y(i))*0.5;
    M0_low(i) = (dM0_low(i+1) + dM0_low(i))*(y(i+1) - y(i))*0.5;
    M_lift(i) = (dM_lift(i+1) + dM_lift(i))*(y(i+1) - y(i))*0.5;
    M_lift_low(i) = (dM_lift_low(i+1) + dM_lift_low(i))*(y(i+1) - y(i))*0.5;
    M_lift_A(i) = (dM_lift_A(i+1) + dM_lift_A(i))*(y(i+1) - y(i))*0.5;
    M_lift_A_low(i) = (dM_lift_A_low(i+1) + dM_lift_A_low(i))*(y(i+1) - y(i))*0.5;
    M_lift_ZF(i) = (dM_lift_ZF(i+1) + dM_lift_ZF(i))*(y(i+1) - y(i))*0.5;
    M_lift_ZF_low(i) = (dM_lift_ZF_low(i+1) + dM_lift_ZF_low(i))*(y(i+1) - y(i))*0.5;
end


M_fuel = zeros(1,length(y));
M_wing = zeros(1,length(y));
M_fuel_low = zeros(1,length(y));
M_wing_low = zeros(1,length(y));
M_fuel_A = zeros(1,length(y));
M_wing_A = zeros(1,length(y));
M_fuel_A_low = zeros(1,length(y));
M_wing_A_low = zeros(1,length(y));
for i = 1:length(y)-1
    M_fuel(i) = (W_fuel(i)-W_fuel(i+1)).*c(i)*(x_FA - x_CG);
    M_wing(i) = (W_wing(i)-W_wing(i+1)).*c(i)*(x_FA - x_CG);
    M_fuel_low(i) = (W_fuel_low(i)-W_fuel_low(i+1)).*c(i)*(x_FA - x_CG);
    M_wing_low(i) = (W_wing_low(i)-W_wing_low(i+1)).*c(i)*(x_FA - x_CG);
    M_fuel_A(i) = (W_fuel_A(i)-W_fuel_A(i+1)).*c(i)*(x_FA - x_CG);
    M_wing_A(i) = (W_wing_A(i)-W_wing_A(i+1)).*c(i)*(x_FA - x_CG);
    M_fuel_A_low(i) = (W_fuel_A_low(i)-W_fuel_A_low(i+1)).*c(i)*(x_FA - x_CG);
    M_wing_A_low(i) = (W_wing_A_low(i)-W_wing_A_low(i+1)).*c(i)*(x_FA - x_CG);
end

M_engine = zeros(1,length(y));
M_engine(length(y)*0.3) = 14601*1.28*n - 2.0193*63000;
M_engine_OEI = zeros(1,length(y));
M_engine_OEI(length(y)*0.3) = 14601*1.28*n;
M_engine_low = zeros(1,length(y));
M_engine_low(length(y)*0.3) = 14601*1.28*n_low - 2.0193*63000; 
M_engine_A = zeros(1,length(y));
M_engine_A(length(y)*0.3) = 14601*1.28*n_A - 2.0193*63000;
M_engine_A_low = zeros(1,length(y));
M_engine_A_low(length(y)*0.3) = 14601*1.28*n_A_low - 2.0193*63000;

dT = M0 + M_lift + M_fuel + M_wing + M_engine;
dT_OEI = M0 + M_lift + M_fuel + M_wing + M_engine_OEI;
dT_A = M0_A + M_lift_A + M_fuel_A + M_wing_A + M_engine_A;
dT_A_low = M0_A_low + M_lift_A_low + M_fuel_A_low + M_wing_A_low + M_engine_A_low;
dT_low = M0_low + M_lift_low + M_fuel_low + M_wing_low + M_engine_low;
dT_ZF = M0 + M_lift_ZF + M_wing + M_engine;
dT_ZF_low = M0_low + M_lift_ZF_low + M_wing_low + M_engine_low;

Torque = zeros(1,length(y));
Torque_OEI = zeros(1,length(y));
Torque_A = zeros(1,length(y));
Torque_A_low = zeros(1,length(y));
Torque_low = zeros(1,length(y));
Torque_ZF = zeros(1,length(y));
Torque_ZF_low = zeros(1,length(y));

for i = 1:length(y)-1
    Torque(i) = sum(dT(i:end));
    Torque_OEI(i) = sum(dT_OEI(i:end));
    Torque_A(i) = sum(dT_A(i:end));
    Torque_A_low(i) = sum(dT_A_low(i:end));
    Torque_low(i) = sum(dT_low(i:end));
    Torque_ZF(i) = sum(dT_ZF(i:end));
    Torque_ZF_low(i) = sum(dT_ZF_low(i:end));
end

Torque_total = Torque*cosd(sweep_FA) + (BM_aero-BM)*sind(sweep_FA);
BM_total = (BM_aero-BM) + Torque*sind(sweep_FA);

Torque_total_OEI = Torque_OEI*cosd(sweep_FA) + (BM_aero-BM)*sind(sweep_FA);
BM_total_OEI = (BM_aero-BM) + Torque_OEI*sind(sweep_FA);

Torque_total_low = Torque_low*cosd(sweep_FA) + (BM_aero_low-BM_low)*sind(sweep_FA);
BM_total_low = (BM_aero_low-BM_low) + Torque_low*sind(sweep_FA);

Torque_total_A = Torque_A*cosd(sweep_FA) + (BM_aero_A-BM_A)*sind(sweep_FA);
BM_total_A = (BM_aero_A-BM_A) + Torque_A*sind(sweep_FA);

Torque_total_A_low = Torque_A_low*cosd(sweep_FA) + (BM_aero_A_low-BM_A_low)*sind(sweep_FA);
BM_total_A_low = (BM_aero_A_low-BM_A_low) + Torque_A_low*sind(sweep_FA);


figure
plot(y,SF_aero - SF,'b-','LineWidth',1.5)
hold on
plot(y,SF_aero_low - SF_low,'b-.','LineWidth',1.5)
plot(y,SF_aero_A - SF_A,'k-','LineWidth',1.5)
plot(y,SF_aero_A_low - SF_A_low,'k-.','LineWidth',1.5)
% plot(y,SF_aero_ZF - SF_ZF,'k-','LineWidth',1.5)
% plot(y,SF_aero_ZF_low - SF_ZF_low,'k-.','LineWidth',1.5)
% plot(y,SF_aero_ZF_A - SF_ZF_A,'r-.','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Shear Force [N]','Interpreter', 'latex', 'FontSize',18)
legend('Load Case 1 ($V_D$ at n = 3.75)','Load Case 1 ($V_D$ at n = -1.5)','Load Case 1 ($V_A$ at n = 2.5)','Load Case 1 ($V_A$ at n = -1)','Interpreter', 'latex', 'FontSize',16)

figure
plot(y,BM_total,'b-','LineWidth',1.5)
hold on
plot(y,BM_total_low,'b-.','LineWidth',1.5)
plot(y,BM_total_A,'k-','LineWidth',1.5)
plot(y,BM_total_A_low,'k-.','LineWidth',1.5)
% plot(y,BM_aero_ZF - BM_ZF,'k-','LineWidth',1.5)
% plot(y,BM_aero_ZF_low - BM_ZF_low,'k-.','LineWidth',1.5)
% plot(y,BM_aero_ZF_A - BM_ZF_A,'r-.','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Bending Moment [Nm]','Interpreter', 'latex', 'FontSize',18)
legend('Load Case 1 ($V_D$ at n = 3.75)','Load Case 1 ($V_D$ at n = -1.5)','Load Case 1 ($V_A$ at n = 2.5)','Load Case 1 ($V_A$ at n = -1)','Interpreter', 'latex', 'FontSize',16)

figure
set(gca,'FontSize',50)
plot(y,Torque_total,'b-','LineWidth',1.5)
hold on
plot(y,Torque_total_low,'b-.','LineWidth',1.5)
plot(y,Torque_total_A,'k-','LineWidth',1.5)
plot(y,Torque_total_A_low,'k-.','LineWidth',1.5)
plot(y,Torque_total_OEI,'g-','LineWidth',1.5)
% plot(y,Torque_ZF,'k-','LineWidth',1.5)
% plot(y,Torque_ZF_low,'k-.','LineWidth',1.5)
hold off
box on 
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Torque [Nm]','Interpreter', 'latex', 'FontSize',18)
legend('Load Case 1 ($V_D$ at n = 3.75)','Load Case 1 ($V_D$ at n = -1.5)','Load Case 1 ($V_A$ at n = 2.5)','Load Case 1 ($V_A$ at n = -1)','Load Case 2 (OEI)','Interpreter', 'latex', 'FontSize',16)



figure
plot(y,dL_dist,'b-','LineWidth',1.5)
hold on
plot(y,W_wing_dist,'r-','LineWidth',1.5)
plot(y,W_fuel_dist,'g-','LineWidth',1.5)
plot(y,W_dist_engine,'k-','LineWidth',1.5)
plot(y,Thrust_dist,'Color',[0.75, 0.75, 0],'LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Load Distribution [N/m]','Interpreter', 'latex', 'FontSize',18)
legend('Lift','Wing Weight','Fuel Weight','Engine Weight','Thrust','Interpreter', 'latex', 'FontSize',16)


%% D-Section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D-Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear y
y = linspace(0,s,100);
c = (c_tip - c_root)*(y/s) + c_root;

SF_dcell = SF_aero(1:100)-SF(1:100);
T_dcell = Torque_total_OEI(1:100);

y_FS_c1 = 0.0989;
y_RS_c1 = 0.0879;

x_FS = 0.15;
x_RS = 0.6;
x_FA = 0.375;

y_FS = y_FS_c1*c;
y_RS = y_RS_c1*c;
y_avg = (y_FS + y_RS)/2;
l_box = c*(x_RS - x_FS);
% 
% E = 73.8*10^9; 
% rho = 2.77*100^3/1000;
% sigma_yield = 430*10^6;
% 
% q_shear = (SF_dcell)./(2*y_avg);
% q_torque = T_dcell./(2.*l_box.*y_avg);
% q_front = abs(q_shear + q_torque);
% q_rear = abs(q_shear - q_torque);
% 
% 
% max_shear_flow = q_front;
% 
% set1_1 = [0,8.1; 1.6,9; 2.7,10; 3.35,11; 3.95,12; 4.4,13; 4.8,14; 5.25,15; 6.05,17; 7.1,20; ...
%        7.4,21; 8,23; 8.4,24; 10.4,31; 11,33; 11.3,34; 11.8,36; 12.3,38; 12.8,40; 14,45; 14.2,46; ...
%        14.9,49; 15.6,52; 15.8,53; 16,54; 16.2,55; 17.1,59];
% p1_1 = polyfit(set1_1(:,1)', set1_1(:,2)', 7);
% 
% % Model curve ratio = 1.5
% set1_15 = [0,6.6; 1.4,7; 2.2,7.5; 2.8,8; 3.55,9; 4.1,10; 4.6,11; 5.1,12; 5.8,13.8; 6.2,15; ...
%        6.8,16.5; 7.3,18; 7.9,20; 8.4,21.5; 9,23.5; 9.2,24; 9.5,25; 9.8,26; 10.8,29.5; 11.2,30.9; ...
%        11.8,33; 12.6,36; 13.6,40; 14.8,45; 15.7,49; 16.8,54; 17.8,59];
% p1_15 = polyfit(set1_15(:,1)', set1_15(:,2)', 7);
% 
% % Model curve ratio = 5
% set1_5 = [0,4.8; 1,5; 2,5.5; 2.6,6; 3.45,7; 4.4,8.5; 5.4,10.5; 6,12; 6.8,14; 8.95,20; 9.8,23; ...
%        10.1,24; 10.4,25; 11.5,29; 12.6,33; 13.6,37; 13.8,38; 14.6,41; 14.8,42; 15.5,45; 16,47; ...
%        17.3,53; 18,56; 18.2,56; 18.4,58; 18.8,60];
% p1_5 = polyfit(set1_5(:,1)', set1_5(:,2)', 7);
% 
% % Plot results
% x1 = 0:0.2:20;
% y1_1 = polyval(p1_1,x1);
% y1_15 = polyval(p1_15,x1);
% y1_5 = polyval(p1_5,x1);
% 
% % figure
% % plot(x1,y1_1,'-b','LineWidth',1.5);
% % hold on
% % plot(x1,y1_15,'-r','LineWidth',1.5);
% % plot(x1,y1_5,'-g','LineWidth',1.5);
% % hold off
% % box on
% % grid on
% % xlim([0 20])
% % ylim([0 60])
% % ax=gca;ax.LineWidth=1;
% % set(gca,'FontSize',15)
% % set(gca,'TickLabelInterpreter','latex')
% % xlabel('$a/sqrt(Rt)$','Interpreter', 'latex', 'FontSize',18)
% % ylabel('$K_s$','Interpreter', 'latex', 'FontSize',18)
% % legend('$1.0$', '$1.5$', '$>5$','Location','southeast','Interpreter', 'latex', 'FontSize',16);
% 
% % Graph 2 (a/b)
% % Model curve ratio = 1
% set2_1 = [0,8.1; 1.6,9; 2.7,10; 3.35,11; 3.95,12; 4.4,13; 4.8,14; 5.25,15; 6.05,17; 7.1,20; ...
%        7.4,21; 8,23; 8.4,24; 10.4,31; 11,33; 11.3,34; 11.8,36; 12.3,38; 12.8,40; 14,45; 14.2,46; ...
%        14.9,49; 15.6,52; 15.8,53; 16,54; 16.2,55; 17.1,59];
% p2_1 = polyfit(set2_1(:,1)', set2_1(:,2)', 7);
% 
% % Model curve ratio = 1.5
% set2_15 = [0,6.6; 1.6,6.8; 2.1,7; 2.8,7.5; 3.35,8; 4,9; 4.55,10; 4.8,10.5; 5.4,11.8; 6.3,14; ...
%        7.1,16; 7.5,17; 8.8,20.5; 10,24; 10.2,24.5; 10.8,26.5; 11,27; 11.8,29.5; 12.4,31.5; 13,33.5; ...
%        14,37; 14.55,39; 15,40.5; 15.7,43; 16.8,47; 17.2,48.5; 17.6,50; 18,51.5; 19.2,56; 19.5,57; ...
%        20,58.8];
% p2_15 = polyfit(set2_15(:,1)', set2_15(:,2)', 7);
% 
% % Model curve ratio = 2
% set2_2 = [0,5.8; 2,6.1; 2.4,6.4; 3,7; 3.4,7.5; 3.85,8; 4.4,9; 4.95,10; 5.4,10.9; 6,12.2; ...
%        6.6,13.5; 7,14.5; 7.6,15.9; 8.45,18; 9.7,21; 10,22; 11,24.5; 11.2,25; 12.2,28; 13.2,31; ...
%        13.8,33; 14.4,35; 15,37; 15.6,39; 16.2,41; 17.2,44.6; 18,47.3; 19,51; 19.6,53; 20,54.5];
% p2_2 = polyfit(set2_2(:,1)', set2_2(:,2)', 7);
% 
% % Model curve ratio = 3
% set2_3 = [0,5.3; 1.4,5.5; 2,5.6; 2.4,5.8; 2.6,6; 3.5,7; 4.6,8.6; 4.8,8.9; 5.4,10; 6,11.2; ...
%        6.5,12; 7.4,14; 7.9,15; 8.8,17; 9.2,18; 10,19.8; 10.5,21; 11.3,23; 12.4,26; 13.2,28; ...
%        14.3,31; 15,33; 15.4,34; 16.4,37; 16.6,37.5; 17.4,39.8; 17.8,41.1; 18.15,42; 18.6,43.5; ...
%        19.1,45; 19.45,46; 19.6,46.5; 19.8,47; 20,47.6];
% p2_3 = polyfit(set2_3(:,1)', set2_3(:,2)', 7);
% 
% % Model curve ratio = inf
% set2_i = [0,4.8; 1.8,5; 2.6,5.5; 3.2,6; 4,7; 4.8,8; 5.4,9; 6,10; 6.6,11; 7.3,12; 7.9,13; ...
%        9.15,15; 9.75,16; 10.4,17; 11,18; 12,19.5; 12.6,20.5; 13.2,21.5; 13.6,22.4; 14.2,23; ...
%        14.8,24; 16,26; 16.4,26.5; 17,27.5; 18,29; 18.55,30; 19.25,31; 19.6,31.5; 19.85,32; 20,32.3];
% p2_i = polyfit(set2_i(:,1)', set2_i(:,2)', 7);
% 
% % Plot results
% x2 = 0:0.2:20;
% y2_1 = polyval(p2_1,x2);
% y2_15 = polyval(p2_15,x2);
% y2_2 = polyval(p2_2,x2);
% y2_3 = polyval(p2_3,x2);
% y2_i = polyval(p2_i,x2);
% 
% % figure
% % plot(x2,y2_1,'-b','LineWidth',1.5);
% % hold on
% % plot(x2,y2_15,'-r','LineWidth',1.5);
% % plot(x2,y2_2,'-g','LineWidth',1.5);
% % plot(x2,y2_3,'-m','LineWidth',1.5);
% % plot(x2,y2_i,'-c','LineWidth',1.5);
% % hold off
% % grid on
% % box on
% % xlim([0 20])
% % ylim([0 60])
% % ax=gca;ax.LineWidth=1;
% % set(gca,'FontSize',15)
% % set(gca,'TickLabelInterpreter','latex')
% % xlabel('$s/sqrt(Rt)$','Interpreter', 'latex', 'FontSize',18)
% % ylabel('$K_s$','Interpreter', 'latex', 'FontSize',18)
% % legend('1.0', '1.5', '2.0', '3.0', 'inf','Location','southeast','Interpreter', 'latex', 'FontSize',16);
% 
% % Save polynomial parameters
% ratiopoly = [p1_1; p1_15; p1_5; p2_1; p2_15; p2_2; p2_3; p2_i];
% 
% tau_tresca = 125*10^6; %%%%%%%%%%%%%%%%%%%%% Update!
% 
% figure
% 
% for k = 1:5 % iterates the whole dcell calculation 3 times to to get converged results
%     n_pseudo_ribs = 0:10;
%     a_pseudo_ribs = s*0.1./((n_pseudo_ribs + 1)*cosd(sweep_LE));
%     t_range = [0.005,0.0055,0.006,0.0065,0.007];
%     t_dcell = t_range(k)*ones(1,length(n_pseudo_ribs));
% 
%     data = [0.14645,  0.03994;0.12408,  0.03795;0.10332,  0.03564;0.08427,  0.03305;0.06699,  0.03023;0.05156,  0.02720;0.03806,  0.02395;0.02653,  0.02039;0.01704,  0.01646;0.00961,  0.01214;0.00428,  0.00767;0.00107,  0.00349;0.0,      0.0];
%     data(:,2) = data(:,2)*1.3;
%     l = 0;
%     for j = 1:length(y)
%         for i = 1:12
%             l(i,j) = sqrt((c(j)*data(i+1,1) - c(j)*data(i,1))^2 + (c(j)*data(i+1,2) - c(j)*data(i,2))^2);
%             A(i,j) = trapz([c(j)*data(i,1),c(j)*data(i+1,1)],[c(j)*data(i,2),c(j)*data(i+1,2)]);
%         end
%         l_curved(j) = (2*sum(l(:,j)) + 2*(0.15-0.14645)*c(j))/2;  
%         A_pseudo_rib(j) = -2*sum(A(:,j));
%     end
% 
%     Radius = (((x_FS.*c).^2)./(y_avg/2)+((y_avg/2).^2)./(x_FS.*c))/2;   % 0.9684
% 
%     for i = 1:length(n_pseudo_ribs)
%         a_b(i,:) = a_pseudo_ribs(i)./l_curved;
%         b_a(i,:) = l_curved/a_pseudo_ribs(i);
%     end
% 
%     for i = 1:length(n_pseudo_ribs)
%         for j = 1:length(y)
%             if a_b(i,j) >= 1
%                 l_sqrt_Rt(i,j) = l_curved(j)/sqrt(Radius(j)*t_dcell(i)); %b_sqrt_Rt
%             else
%                 l_sqrt_Rt(i,j) = a_pseudo_ribs(i)/sqrt(Radius(j)*t_dcell(i)); %a_sqrt_Rt
%             end
%         end
%     end
%     for i = 1:length(n_pseudo_ribs)
%         Ks_dcell(i,:) = DSectionRatioFunc(b_a(i,:),l_sqrt_Rt(i,:),ratiopoly);
%         tau_crit(i,:) = Ks_dcell(i,:).*E.*t_dcell(i).^3./(l_curved.^2);
%         tau_crit_min(i) = min(tau_crit(i,:));
% %         t_panel(i,:) = sqrt(tau_tresca.*l_curved.^2./(Ks_dcell(i,:).*E));
% %         t_opt(i) = max(t_panel(i,:));
%     end
% 
%     plot(n_pseudo_ribs,tau_crit_min,'LineWidth',1.5)
%     hold on
% end
% plot([0 10],[max_shear_flow(1),max_shear_flow(1)],'k--','LineWidth',1.5)
% hold off
% grid on
% box on
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex')
% xlabel('Number of Pseudo Ribs','Interpreter', 'latex', 'FontSize',18)
% ylabel('Minimum Buckling Shear Flow [N/m]','Interpreter', 'latex', 'FontSize',18)
% legend('$t = 5\,mm$','$t = 5.5\,mm$','$t = 6\,mm$','$t = 6.5\,mm$','$t = 7\,mm$','Max Applied Shear Flow','Interpreter', 'latex', 'FontSize',18)
% 
% m_dcell = rho*t_range*a_pseudo_ribs(1)*(l_curved(1)+l_curved(end))/2;
% t_pseudo_rib = 0.001;
% 
% for i = 1:length(n_pseudo_ribs)
%     m_pseudo_ribs(i) = 0;
%     if n_pseudo_ribs(i) > 0
%         for j = 1:n_pseudo_ribs(i)
%             m_pseudo_ribs(i) = m_pseudo_ribs(i) + rho*A_pseudo_rib(floor((length(y)*a_pseudo_ribs(i)*j)/(0.1*s/cosd(sweep_LE))))*t_pseudo_rib*n_pseudo_ribs(i);
%         end
%     end
% end
% 
% for i = 1:length(t_range)
%     for j = 1:length(n_pseudo_ribs)
%         m_total_dcell(i,j) = m_dcell(i) + m_pseudo_ribs(j);
%     end
% end
% 
% figure
% plot(n_pseudo_ribs,m_total_dcell(1,:),'LineWidth',1.5)
% hold on
% plot(n_pseudo_ribs,m_total_dcell(2,:),'LineWidth',1.5)
% plot(n_pseudo_ribs,m_total_dcell(3,:),'LineWidth',1.5)
% plot(n_pseudo_ribs,m_total_dcell(4,:),'LineWidth',1.5)
% plot(n_pseudo_ribs,m_total_dcell(5,:),'LineWidth',1.5)
% hold off
% grid on
% box on
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex')
% xlabel('Number of Pseudo Ribs','Interpreter', 'latex', 'FontSize',18)
% ylabel('Overall D-Cell Mass [kg]','Interpreter', 'latex', 'FontSize',18)
% legend('$t = 5\,mm$','$t = 5.5\,mm$','$t = 6\,mm$','$t = 6.5\,mm$','$t = 7\,mm$','Interpreter', 'latex', 'FontSize',18)


%% Wing Configuration

%%%%%%%%%%%%%%%%%%%%%%%%% Plot Wing Configuration %%%%%%%%%%%%%%%%%%%%%%%%

load ribLoc2.mat

n = 16; % number of stringers
n_pseudo = 0;
b_stringers = 0.12081;
figure
p1 = plot([0,s,s,0,0],[c_root,c_root-s*tand(sweep_LE),c_root-s*tand(sweep_LE)-c_tip,0,c_root],'k-','LineWidth',1.5);
hold on
z_LE = @(x) ((c_root-(c_root-s*tand(sweep_LE)))/(0-s))*x + c_root;
z = @(x) (((c_root*(1-x_FS)-(c_tip*(1-x_FS)+c_root-s*tand(sweep_LE)-c_tip))/(0-s))*x + c_root*(1-x_FS))*(x>=1.5) + (((c_root*(1-x_FS)-(c_tip*(1-x_FS)+c_root-s*tand(sweep_LE)-c_tip))/(0-s)-0.15)*x + (c_root*(1-x_FS)+0.215))*(x<1.5);
z_RS = @(x) (((c_root*(1-x_RS)-(c_tip*(1-x_RS)+c_root-s*tand(sweep_LE)-c_tip))/(0-s))*x + c_root*(1-x_RS))*(x>=1.5) + (((c_root*(1-x_RS)-(c_tip*(1-x_RS)+c_root-s*tand(sweep_LE)-c_tip))/(0-s)+0.15)*x + (c_root*(1-x_RS)-0.215))*(x<1.5);
for i = 1:n
    hoe = c_tip*(1-x_RS)+c_root-s*tand(sweep_LE)-c_tip+i*b_stringers;
    if hoe <= c_tip*(1-x_FS)+c_root-s*tand(sweep_LE)-c_tip
        x_st(i,:) = [0,s];
        y_st(i,:) = [c_root*(1-x_RS)+i*b_stringers,hoe];
        p3 = plot(x_st(i,:),y_st(i,:),'color',[0, 0.40, 0],'LineWidth',1.5);
    else
        z_stringer = @(x) ((c_root*(1-x_RS)+i*b_stringers-hoe)/(0-s))*x + c_root*(1-x_RS)+i*b_stringers;
        Int = fzero(@(x) z(x)-z_stringer(x), 1);
        plot([0,Int],[c_root*(1-x_RS)+i*b_stringers,z(Int)],'color',[0, 0.40, 0],'LineWidth',1.5)
    end
    plot([0,s],[c_root*(1-x_RS),c_tip*(1-x_RS)+c_root-s*tand(sweep_LE)-c_tip],'color',[0, 0.40, 0],'LineWidth',1.5)
    plot([0,0.75],[c_root*(1-x_RS-0.03),z_RS(0.75)],'color',[0, 0.40, 0],'LineWidth',1.5)
end

for j = 2:length(ribLoc2)-1
    p4 = plot([0+ribLoc2(j),0+ribLoc2(j)],[z(ribLoc2(j)),z_RS(ribLoc2(j))],'r-','LineWidth',1.5); %z(ribLoc2(j))-l_box(round(ribLoc2(j)*100/s)
end

if n_pseudo > 0
    for k = 1:n_pseudo
        spacing = floor(length(y)/(n_pseudo+1));
        p5 = plot([0+y(spacing*k),0+y(spacing*k)],[z_LE(y(spacing*k)),z(y(spacing*k))],'Color',[0.75, 0.75, 0],'LineWidth',1.5);
    end
end
for i = 1:length(y)-1
    p2 = plot([y(i),y(i+1)],[z(y(i)),z(y(i+1))],'b-','LineWidth',1.5)%plot([0,s],[c_root*(1-x_FS),c_tip*(1-x_FS)+c_root-s*tand(sweep_LE)-c_tip],'b-','LineWidth',1.5);
    p2 = plot([y(i),y(i+1)],[z_RS(y(i)),z_RS(y(i+1))],'b-','LineWidth',1.5)
end
p6 = fill([s*0.1,s*0.7,s*0.7,s*0.1,s*0.1],[c_root*0.25-s*0.1*tand(sweep_hinge),c_root*0.25-s*0.7*tand(sweep_hinge),c_root*0.25-s*0.7*tand(sweep_hinge)-0.57,c_root*0.25-s*0.1*tand(sweep_hinge)-0.93,c_root*0.25-s*0.1*tand(sweep_hinge)],[0, 0.30, 0.40]); %0 0.9 0.5
p7 = fill([s*0.7,s*0.95,s*0.95,s*0.7,s*0.7],[c_root*0.25-s*0.7*tand(sweep_hinge),c_root*0.25-s*0.95*tand(sweep_hinge),c_root*0.25-s*0.95*tand(sweep_hinge)-0.44,c_root*0.25-s*0.7*tand(sweep_hinge)-0.57,c_root*0.25-s*0.7*tand(sweep_hinge)],'k'); %0 0.9 0.5
p8 = fill([s*0.1,s*0.99,s*0.99,s*0.1,s*0.1],[c_root-s*0.1*tand(sweep_LE),c_root-s*0.99*tand(sweep_LE),c_root*0.9-s*0.99*tand(sweep_10),c_root*0.9-s*0.1*tand(sweep_10),c_root-s*0.1*tand(sweep_LE)],[0.827,0.827,0.827]);
hold off
box on 
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise Position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Chordwise Position [m]','Interpreter', 'latex', 'FontSize',18)
legend([p2 p3 p4 p8 p6 p7],{'Spars','Stringers','Ribs','Slats','Flaps','Ailerons'},'Interpreter', 'latex', 'FontSize',16)


