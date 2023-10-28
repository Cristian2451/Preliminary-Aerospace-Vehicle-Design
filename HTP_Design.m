clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_FA = 0.375; % Flexural Axis 
x_AC = 0.25; % Aerodynamic Centre
x_CG = 0.4; % Wing + Fuel Centre of Gravity (update to more accurate one!!)

MTOW = 385434.9; % [N]
n = 3.75;
n_low = -1.5;
n_A = 2.5;
n_A_low = -1;
Sref_h = 10;
Sref_w = 79.471;
taper = 0.45;
AR = 4;
sweep = 17;
rho = 1.225;
sweep_LE = atand(tand(sweep)-(4/AR)*((0 - x_AC)*(1-taper)/(1+taper))); % at quarter chord
sweep_FA = atand(tand(sweep)-(4/AR)*((x_FA - x_AC)*(1-taper)/(1+taper)));
sweep_hinge = atand(tand(sweep)-(4/AR)*((0.75 - x_AC)*(1-taper)/(1+taper)));

% Analysing HTP from fuselage FoR

b = sqrt(AR*Sref_h);
s = b*0.5; %/cosd(sweep); % Update sweep to flexural axis!!!!!!!
y = linspace(0,s,1000);
c_root = 2*Sref_h/(b*(1+taper));
c_tip = taper*c_root; 

c_root_w = 4.0259;

W_h = 0.5*155.94*9.81*n;
W_h_low = 0.5*155.94*9.81*n_low;
W_h_land = 0.5*155.94*9.81;
W_h_A = 0.5*155.94*9.81*n_A;
W_h_A_low = 0.5*155.94*9.81*n_A_low;

W_h_root = (2*W_h)/(s*(c_tip/c_root + 1)); %Weight at root (N/m)
W_h_tip = (c_tip/c_root)*W_h_root; %Weight at tip (N/m)
W_dist = ((W_h_tip-W_h_root)/s).*y + W_h_root; %Weight distribution along span

W_h_root_A = (2*W_h_A)/(s*(c_tip/c_root + 1)); %Weight at root (N/m)
W_h_tip_A = (c_tip/c_root)*W_h_root_A; %Weight at tip (N/m)
W_dist_A = ((W_h_tip_A-W_h_root_A)/s).*y + W_h_root_A; %Weight distribution along span

W_h_root_A_low = (2*W_h_A_low)/(s*(c_tip/c_root + 1)); %Weight at root (N/m)
W_h_tip_A_low = (c_tip/c_root)*W_h_root_A_low; %Weight at tip (N/m)
W_dist_A_low = ((W_h_tip_A_low-W_h_root_A_low)/s).*y + W_h_root_A_low; %Weight distribution along span

W_h_root_low = (2*W_h_low)/(s*(c_tip/c_root + 1)); %Weight at root (N/m)
W_h_tip_low = (c_tip/c_root)*W_h_root_low; %Weight at tip (N/m)
W_dist_low = ((W_h_tip_low-W_h_root_low)/s).*y + W_h_root_low; %Weight distribution along span

W_h_root_land = (2*W_h_land)/(s*(c_tip/c_root + 1)); %Weight at root (N/m)
W_h_tip_land = (c_tip/c_root)*W_h_root_land; %Weight at tip (N/m)
W_dist_land = ((W_h_tip_land-W_h_root_land)/s).*y + W_h_root_land; %Weight distribution along span

%%%%%%%%%%%%%%%%%%%%%%%%%%% Inertia Loads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = zeros(1,length(y));
W_A = zeros(1,length(y));
W_A_low = zeros(1,length(y));
W_low = zeros(1,length(y));
W_land = zeros(1,length(y));

for i = 1:length(y)-1
    W(i) = (W_dist(i+1) + W_dist(i))*(y(i+1) - y(i))*0.5;
    W_A(i) = (W_dist_A(i+1) + W_dist_A(i))*(y(i+1) - y(i))*0.5;
    W_A_low(i) = (W_dist_A_low(i+1) + W_dist_A_low(i))*(y(i+1) - y(i))*0.5;
    W_low(i) = (W_dist_low(i+1) + W_dist_low(i))*(y(i+1) - y(i))*0.5;
    W_land(i) = (W_dist_land(i+1) + W_dist_land(i))*(y(i+1) - y(i))*0.5;
end

SF = zeros(1,length(y));
dM = zeros(1,length(y));
BM = zeros(1,length(y));
SF_A = zeros(1,length(y));
dM_A = zeros(1,length(y));
BM_A = zeros(1,length(y));
SF_A_low = zeros(1,length(y));
dM_A_low = zeros(1,length(y));
BM_A_low = zeros(1,length(y));
SF_low = zeros(1,length(y));
dM_low = zeros(1,length(y));
BM_low = zeros(1,length(y));
SF_land = zeros(1,length(y));
dM_land = zeros(1,length(y));
BM_land = zeros(1,length(y));

for i = 1:length(y)-1
    SF(i) = sum(W(i:end));
    SF_A(i) = sum(W_A(i:end));
    SF_A_low(i) = sum(W_A_low(i:end));
    SF_low(i) = sum(W_low(i));
    SF_land(i) = sum(W_land(i));
end

for i = 1:length(y)-1
    dM(i) = (SF(i+1) + SF(i))*(y(i+1) - y(i))*0.5;
    dM_A(i) = (SF_A(i+1) + SF_A(i))*(y(i+1) - y(i))*0.5;
    dM_A_low(i) = (SF_A_low(i+1) + SF_A_low(i))*(y(i+1) - y(i))*0.5;
    dM_low(i) = (SF_low(i+1) + SF_low(i))*(y(i+1) - y(i))*0.5;
    dM_land(i) = (SF_land(i+1) + SF_land(i))*(y(i+1) - y(i))*0.5;
end

for i = 1:length(y)-1
    BM(i) = sum(dM(i:end));
    BM_A(i) = sum(dM_A(i:end));
    BM_A_low(i) = sum(dM_A_low(i:end));
    BM_low(i) = sum(dM_low(i:end));
    BM_land(i) = sum(dM_land(i:end));
end

% figure
% plot(y,SF,'b-','LineWidth',1.5)
% grid on
% box on
% xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',14)
% ylabel('Inertia Shear Force [N]','Interpreter', 'latex', 'FontSize',14)
% legend('Load Case 1 (MTOW)','Interpreter', 'latex', 'FontSize',12)
% 
% figure
% plot(y,BM,'b-','LineWidth',1.5)
% grid on
% box on
% xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',14)
% ylabel('Inertia Bending Moment [Nm]','Interpreter', 'latex', 'FontSize',14)
% legend('Load Case 1 (MTOW)','Interpreter', 'latex', 'FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Aero Loads %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_h = 28.739 + 0.25*c_root;
x_w = 15.14;
x_cg = 14.84;

Cm0 = -0.023;
V_ne = 196.884;
V_A = 147.5;

M0 = 0.5*rho*V_ne^2*Sref_w*2.991*Cm0;
M0_A = 0.5*rho*V_A^2*Sref_w*2.991*Cm0;


% Find lift at horizontal tailplane by taking moments from x_cg
L_w = n*MTOW;
L_h = (L_w*(x_cg - x_w) + M0)/(x_h - x_cg);
L0_h = 2*L_h/(pi*s);
dL_h = L0_h*sqrt(1-(y/s).^2);

L_w_A = n_A*MTOW;
L_h_A = (L_w_A*(x_cg - x_w) + M0_A)/(x_h - x_cg);
L0_h_A = 2*L_h_A/(pi*s);
dL_h_A = L0_h_A*sqrt(1-(y/s).^2);

L_w_A_low = n_A_low*MTOW;
L_h_A_low = (L_w_A_low*(x_cg - x_w))/(x_h - x_cg);
L0_h_A_low = 2*L_h_A_low/(pi*s);
dL_h_A_low = L0_h_A_low*sqrt(1-(y/s).^2);

L_w_low = n_low*MTOW;
L_h_low = (L_w_low*(x_cg - x_w))/(x_h - x_cg);
L0_h_low = 2*L_h_low/(pi*s);
dL_h_low = L0_h_low*sqrt(1-(y/s).^2);

x_h = 28.739 + 0.25*c_root;
x_cg = 15.08;
x_mlg = 17.71;

% F_h + F_mlg = MTOW (sum of vertical forces = 0)

% F_h*(x_h - x_cg) + F_mlg*(x_mlg - x_cg) = 0 (sum of moments = 0)

% Solving simultaneous equation gives:

F_h = -MTOW*(x_mlg - x_cg)/((x_h - x_cg) - (x_mlg - x_cg));
F_mlg = MTOW - F_h;

L_h_land = F_h; %-58413;
L0_h_land = 2*L_h_land/(pi*s);
dL_h_land = L0_h_land*sqrt(1-(y/s).^2);

Lift = zeros(1,length(y));
SF_aero = zeros(1,length(y));
dM_aero = zeros(1,length(y));
BM_aero = zeros(1,length(y));

Lift_A = zeros(1,length(y));
SF_aero_A = zeros(1,length(y));
dM_aero_A = zeros(1,length(y));
BM_aero_A = zeros(1,length(y));

Lift_A_low = zeros(1,length(y));
SF_aero_A_low = zeros(1,length(y));
dM_aero_A_low = zeros(1,length(y));
BM_aero_A_low = zeros(1,length(y));

Lift_low = zeros(1,length(y));
SF_aero_low = zeros(1,length(y));
dM_aero_low = zeros(1,length(y));
BM_aero_low = zeros(1,length(y));

Lift_land = zeros(1,length(y));
SF_aero_land = zeros(1,length(y));
dM_aero_land = zeros(1,length(y));
BM_aero_land = zeros(1,length(y));

for i = 1:length(y)-1
    Lift(i) = (dL_h(i+1) + dL_h(i))*(y(i+1) - y(i))*0.5;
    Lift_A(i) = (dL_h_A(i+1) + dL_h_A(i))*(y(i+1) - y(i))*0.5;
    Lift_A_low(i) = (dL_h_A_low(i+1) + dL_h_A_low(i))*(y(i+1) - y(i))*0.5;
    Lift_low(i) = (dL_h_low(i+1) + dL_h_low(i))*(y(i+1) - y(i))*0.5;
    Lift_land(i) = (dL_h_land(i+1) + dL_h_land(i))*(y(i+1) - y(i))*0.5;
end

for i = 1:length(y)-1
    SF_aero(i) = sum(Lift(i:end));
    SF_aero_A(i) = sum(Lift_A(i:end));
    SF_aero_A_low(i) = sum(Lift_A_low(i:end));
    SF_aero_low(i) = sum(Lift_low(i:end));
    SF_aero_land(i) = sum(Lift_land(i:end));
end

for i = 1:length(y)-1
    dM_aero(i) = (SF_aero(i+1) + SF_aero(i))*(y(i+1) - y(i))*0.5;
    dM_aero_A(i) = (SF_aero_A(i+1) + SF_aero_A(i))*(y(i+1) - y(i))*0.5;
    dM_aero_A_low(i) = (SF_aero_A_low(i+1) + SF_aero_A_low(i))*(y(i+1) - y(i))*0.5;
    dM_aero_low(i) = (SF_aero_low(i+1) + SF_aero_low(i))*(y(i+1) - y(i))*0.5;
    dM_aero_land(i) = (SF_aero_land(i+1) + SF_aero_land(i))*(y(i+1) - y(i))*0.5;
end

for i = 1:length(y)-1
    BM_aero(i) = sum(dM_aero(i:end));
    BM_aero_A(i) = sum(dM_aero_A(i:end));
    BM_aero_A_low(i) = sum(dM_aero_A_low(i:end));
    BM_aero_low(i) = sum(dM_aero_low(i:end));
    BM_aero_land(i) = sum(dM_aero_land(i:end));
end


figure
plot(y,dL_h,'b-','LineWidth',1.5)
hold on
plot(y,-W_dist,'r-','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Load Distribution [N/m]','Interpreter', 'latex', 'FontSize',18)
legend('Lift','Weight','Interpreter', 'latex', 'FontSize',16)


% figure
% plot(y,SF_aero,'b-','LineWidth',1.5)
% grid on
% box on
% xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',14)
% ylabel('Aerodynamic Shear Force [N]','Interpreter', 'latex', 'FontSize',14)
% legend('Load Case 1 (MTOW)','Interpreter', 'latex', 'FontSize',12)
% 
% figure
% plot(y,BM_aero,'b-','LineWidth',1.5)
% grid on
% box on
% xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',14)
% ylabel('Aerodynamic Bending Moment [Nm]','Interpreter', 'latex', 'FontSize',14)
% legend('Load Case 1 (MTOW)','Interpreter', 'latex', 'FontSize',12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Torque %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Moment defined as +ve in anti-clockwise direction

c = (c_tip - c_root)*(y/s) + c_root;

dM_lift = -dL_h.*c*(x_FA - x_AC);
dM_lift_A = -dL_h_A.*c*(x_FA - x_AC);
dM_lift_A_low = -dL_h_A_low.*c*(x_FA - x_AC);
dM_lift_low = -dL_h_low.*c*(x_FA - x_AC);
dM_lift_land = -dL_h_land.*c*(x_FA - x_AC);

M_lift = zeros(1,length(y));
M_lift_A = zeros(1,length(y));
M_lift_A_low = zeros(1,length(y));
M_lift_low = zeros(1,length(y));
M_lift_land = zeros(1,length(y));

for i = 1:length(y)-1
    M_lift(i) = (dM_lift(i+1) + dM_lift(i))*(y(i+1) - y(i))*0.5;
    M_lift_A(i) = (dM_lift_A(i+1) + dM_lift_A(i))*(y(i+1) - y(i))*0.5;
    M_lift_A_low(i) = (dM_lift_A_low(i+1) + dM_lift_A_low(i))*(y(i+1) - y(i))*0.5;
    M_lift_low(i) = (dM_lift_low(i+1) + dM_lift_low(i))*(y(i+1) - y(i))*0.5;
    M_lift_land(i) = (dM_lift_land(i+1) + dM_lift_land(i))*(y(i+1) - y(i))*0.5;
end

M_wing = W.*c*(x_FA - x_CG);
M_wing_A = W_A.*c*(x_FA - x_CG);
M_wing_A_low = W_A_low.*c*(x_FA - x_CG);
M_wing_low = W_low.*c*(x_FA - x_CG);
M_wing_land = W_land.*c*(x_FA - x_CG);

dT = M_lift + M_wing;
dT_A = M_lift_A + M_wing_A;
dT_A_low = M_lift_A_low + M_wing_A_low;
dT_low = M_lift_low + M_wing_low;
dT_land = M_lift_land + M_wing_land;

Torque = zeros(1,length(y));
Torque_A = zeros(1,length(y));
Torque_A_low = zeros(1,length(y));
Torque_low = zeros(1,length(y));
Torque_land = zeros(1,length(y));

for i = 1:length(y)-1
    Torque(i) = sum(dT(i:end));
    Torque_A(i) = sum(dT_A(i:end));
    Torque_A_low(i) = sum(dT_A_low(i:end));
    Torque_low(i) = sum(dT_low(i:end));
    Torque_land(i) = sum(dT_land(i:end));
end

% figure
% plot(y,M0)
% hold on 
% plot(y,M_lift)
% plot(y,M_fuel)
% plot(y,M_wing)
% plot(y,M_engine)
% hold off
% legend('M0','M_lift','M_fuel','M_wing','M_engine')

SF_total = SF_aero - SF;
SF_total_A = SF_aero_A - SF_A;
SF_total_A_low = SF_aero_A_low - SF_A_low;
SF_total_low = SF_aero_low - SF_low;
SF_total_land = SF_aero_land - SF_land;

Torque_total = Torque + (BM_aero - BM)*sind(sweep_FA);
Torque_total_A = Torque_A + (BM_aero_A - BM_A)*sind(sweep_FA);
Torque_total_A_low = Torque_A_low + (BM_aero_A_low - BM_A_low)*sind(sweep_FA);
Torque_total_low = Torque_low + (BM_aero_low - BM_low)*sind(sweep_FA);
Torque_total_land = Torque_land + (BM_aero_land - BM_land)*sind(sweep_FA);

BM_tot = (BM_aero - BM) + Torque*sind(sweep_FA); %+ Torque*cosd(sweep_FA);
BM_tot_A = (BM_aero_A - BM_A) + Torque_A*sind(sweep_FA);
BM_tot_A_low = (BM_aero_A_low - BM_A_low) + Torque_A_low*sind(sweep_FA);
BM_total_low = (BM_aero_low - BM_low) + Torque_low*sind(sweep_FA); %+ Torque_low*cosd(sweep_FA);   % check if it is cos or sin
BM_total_land = (BM_aero_land - BM_land) + Torque_land*sind(sweep_FA); %+ Torque_land*cosd(sweep_FA);

figure
plot(y,SF_total,'b-','LineWidth',1.5)
hold on
plot(y,SF_total_low,'b-.','LineWidth',1.5)
plot(y,SF_total_A,'k-','LineWidth',1.5)
plot(y,SF_total_A_low,'k-.','LineWidth',1.5)
plot(y,SF_total_land,'r-','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Shear Force [N]','Interpreter', 'latex', 'FontSize',18)
legend('Load Case 1 ($V_D$ at n = 3.75)','Load Case 1 ($V_D$ at n = -1.5)','Load Case 1 ($V_A$ at n = 2.5)','Load Case 1 ($V_A$ at n = -1)','Load Case 3','Interpreter', 'latex', 'FontSize',16)

figure
plot(y,BM_tot,'b-','LineWidth',1.5)
hold on
plot(y,BM_total_low,'b-.','LineWidth',1.5)
plot(y,BM_tot_A,'k-','LineWidth',1.5)
plot(y,BM_tot_A_low,'k-.','LineWidth',1.5)
plot(y,BM_total_land,'r-','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Bending Moment [Nm]','Interpreter', 'latex', 'FontSize',18)
legend('Load Case 1 ($V_D$ at n = 3.75)','Load Case 1 ($V_D$ at n = -1.5)','Load Case 1 ($V_A$ at n = 2.5)','Load Case 1 ($V_A$ at n = -1)','Load Case 3','Interpreter', 'latex', 'FontSize',16)

figure
plot(y,Torque_total,'b-','LineWidth',1.5)
hold on
plot(y,Torque_total_low,'b-.','LineWidth',1.5)
plot(y,Torque_total_A,'k-','LineWidth',1.5)
plot(y,Torque_total_A_low,'k-.','LineWidth',1.5)
plot(y,Torque_total_land,'r-','LineWidth',1.5)
hold off
box on 
grid on
%ylim([-6000 0])
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Torque [Nm]','Interpreter', 'latex', 'FontSize',18)
legend('Load Case 1 ($V_D$ at n = 3.75)','Load Case 1 ($V_D$ at n = -1.5)','Load Case 1 ($V_A$ at n = 2.5)','Load Case 1 ($V_A$ at n = -1)','Load Case 3','Interpreter', 'latex', 'FontSize',16)


%% Design 

% Parameters

E = 73.8*10^9; % 2024-T351
rho = 2.77*100^3/1000;
sigma_yield = 430*10^6;

x_FS = 0.15;
x_RS = 0.6;
y_FS_c1 = 0.08038; % From airfoil tools (For c = 1)
y_RS_c1 = 0.06853; % From airfoil tools (For c = 1)
y_avg_c1 = (y_RS_c1+y_FS_c1)/2;

y_FS = y_FS_c1*c;
y_RS = y_RS_c1*c;
y_avg = (y_FS + y_RS)/2;
l_box = c*(x_RS - x_FS);

x_airfoil = [1.00, 0.99572, 0.98296, 0.96194, 0.93301, 0.89668, 0.85355, 0.80438, 0.75, 0.69134, 0.62941, 0.56526, 0.5, 0.43474, 0.37059, 0.33928, 0.30866, 0.27886, 0.25, 0.22221, 0.19562, 0.17033, 0.14645, 0.12408, 0.10332, 0.08427, 0.06699, 0.05156, 0.03806, 0.02653, 0.01704, 0.00961, 0.00428, 0.00107, 0, 0.00107, 0.00428, 0.00961, 0.01704, 0.02653, 0.03806, 0.05156, 0.06699, 0.08427, 0.10332, 0.12408, 0.14645, 0.17033, 0.19562, 0.22221, 0.25, 0.27886, 0.30866, 0.33928, 0.37059, 0.43474, 0.5, 0.56526, 0.62941, 0.69134, 0.75, 0.80438, 0.85355, 0.89668, 0.93301, 0.96194, 0.98296, 0.99572, 1];
y_airfoil = [0.0, 0.00057, 0.00218, 0.00463, 0.0077, 0.01127, 0.01522, 0.01945, 0.02384, 0.02823, 0.03247, 0.03638, 0.03978, 0.04248, 0.04431, 0.04484, 0.04509, 0.04504, 0.04466, 0.04397, 0.04295, 0.04161, 0.03994, 0.03795, 0.03564, 0.03305, 0.03023, 0.0272, 0.02395, 0.02039, 0.01646, 0.01214, 0.00767, 0.00349, 0, -0.00349, -0.00767, -0.01214, -0.01646, -0.02039, -0.02395, -0.0272, -0.03023, -0.03305, -0.03564, -0.03795, -0.03994, -0.04161, -0.04295, -0.04397, -0.04466, -0.04504, -0.04509, -0.04484, -0.04431, -0.04248, -0.03978, -0.03638, -0.03247, -0.02823, -0.02384, -0.01945, -0.01522, -0.01127, -0.0077, -0.00463, -0.00218, -0.00057, 0];

figure 
plot([x_FS,x_FS,x_RS,x_RS,x_FS],[y_avg_c1/2,-y_avg_c1/2,-y_avg_c1/2,y_avg_c1/2,y_avg_c1/2],'b-','LineWidth',2)
hold on
plot(x_FA,0,'b*','LineWidth',2,'MarkerSize',14)
%plot(x_CG,0,'ko','LineWidth',2,'MarkerSize',14)
plot(x_airfoil,y_airfoil,'k-','LineWidth',2)
hold off
grid minor
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
xlabel('x/c [m]','Interpreter', 'latex', 'FontSize',24)
ylabel('y/c [m]','Interpreter', 'latex', 'FontSize',24)
legend('Wingbox','Flexural Axis','Interpreter', 'latex', 'FontSize',22)

BM_total = -BM_total_land; %-(BM_aero_land - BM_land);

N_root = abs((BM_total(1))./(l_box(1).*y_avg(1))); %compressive load per unit length

%%%%%%%%%%%%%%%%%%%%% stringers (t_skin constant) %%%%%%%%%%%%%%%%%%%%%%%%%

n_stringers = 0:40;

b_stringers = l_box(1)./(n_stringers + 1);

% Values to optimize (for minimum mass): n° of stringers, t_s, h
% Keeping sigma_crit = N/t

d_h_CT = 0.3; % Value indicated on ESDU plot

ts_t_CT = 1; % selected for max farrar  efficiency
%As_bt_CT = 1.4; % selected for max farrar efficiency


% If As_bt_CT is kept constant, the area and hence height of stringers will
% vary with number of stringers.
% It's more reasonable to set stringer height (30% of winbox height at tip)
% hence As_bt_CT will vary with number of stringers and so will do Farrar.

%b_h_CT = (1+2*d_h_CT)*ts_t_CT/As_bt_CT; 
%h_b_CT = 1/b_h_CT;
h_CT = round(0.4*y_RS(end)*1000)/1000;  % 40% of rear spar height at tip
h_b_CT = h_CT./b_stringers;

hb_1 = [0.101663586, 0.12754159, 0.153419593, 0.22181146, 0.273567468, 0.336414048, 0.410351201, 0.46025878, 0.547134935, 0.600739372, 0.641404806, 0.689463956, 0.748613678, 0.802218115, 0.859519409, 0.896487985, 0.926062847, 0.959334566, 1];
K_1 = [5.488946459, 5.383592401, 5.278238342, 5.056994819, 4.920034542, 4.814680484, 4.709326425, 4.677720207, 4.603972366, 4.561830743, 4.519689119, 4.46701209, 4.393264249, 4.308981002, 4.214162349, 4.129879102, 4.066666667, 3.961312608, 3.813816926];

[xData4, yData4] = prepareCurveData(hb_1, K_1);
f1 = fit(xData4, yData4,'pchip');

t_skin_CT = zeros(1,length(n_stringers));

for i = 1:length(n_stringers)
    t_skin_CT(i) = ((N_root.*b_stringers(i).^2)./(f1(h_b_CT(i))*E)).^(1/3);
end

t_skin_CT(t_skin_CT < 0.001) = 0.001;
%t_skin_CT = ceil(t_skin_CT*1000)/1000;

% A_stringer_CT = As_bt_CT.*b_stringers.*t_skin_CT;
A_stringer_CT = h_CT*(1 + 2*d_h_CT)*ts_t_CT*t_skin_CT;
As_bt_CT = A_stringer_CT./(b_stringers.*t_skin_CT);

figure
plot(n_stringers,As_bt_CT)
xlabel('Number of stringers','Interpreter', 'latex', 'FontSize',18)
ylabel('As_bt_CT','Interpreter', 'latex', 'FontSize',18)

t_eff_CT = t_skin_CT + A_stringer_CT./b_stringers;

for i = 1:length(y)
    A_eff_CT(:,i) = t_skin_CT.*l_box(i) + A_stringer_CT.*(floor(l_box(i)./b_stringers));
end

for i = 1:length(n_stringers)
    V_panel_CT(i) = 0;
    for j = 1:length(y)-1
        V_panel_CT(i) = V_panel_CT(i) + (A_eff_CT(i,j) + A_eff_CT(i,j+1))*0.5*y(2);
    end
end

m_panel_CT = V_panel_CT*rho;

%%%%%%%%%%%%%%%%%%%%% stringers (t_skin varying) %%%%%%%%%%%%%%%%%%%%%%%%%

d_h = 0.3; % Value indicated on ESDU plot
ts_t = 0.5;
As_bt = 0.6;

h = round(0.4*y_RS(end)*1000)/1000; 
h_b = h./b_stringers;

hb_0_5 = [0.099815157, 0.134935305, 0.171903882, 0.208872458, 0.258780037, 0.323475046, 0.384473198, 0.467652495, 0.537892791, 0.558225508, 0.569316081, 0.589648799, 0.609981516, 0.630314233, 0.654343808, 0.676524954, 0.704251386, 0.731977819, 0.763401109, 0.789279113, 0.831792976, 0.870609982, 0.918669131, 0.972273567, 0.998151571];
K_0_5 = [4.150949914, 4.056131261, 3.98238342, 3.940241796, 3.898100173, 3.855958549, 3.80328152, 3.740069085, 3.666321244, 3.645250432, 3.529360967, 3.360794473, 3.181692573, 3.013126079, 2.823488774, 2.65492228, 2.454749568, 2.275647668, 2.107081174, 1.959585492, 1.769948187, 1.611917098, 1.443350604, 1.306390328, 1.232642487];

[xData9, yData9] = prepareCurveData(hb_0_5, K_0_5);
f0_5 = fit(xData9, yData9,'pchip');

N = abs((BM_total)./(l_box.*y_avg)); %compressive load per unit length

for i = 1:length(n_stringers)
    t_skin(i,:) = ((N.*b_stringers(i).^2)./(f0_5(h_b(i))*E)).^(1/3);
end

t_skin(t_skin < 0.001) = 0.001;
t_skin_round = ceil(t_skin*1000)/1000;

figure
plot(y,t_skin(1,:),'b-','LineWidth',1.5)
hold on
plot(y,t_skin(10+1,:),'r-','LineWidth',1.5)
plot(y,t_skin(25+1,:),'g-','LineWidth',1.5)
plot(y,t_skin(40+1,:),'k-','LineWidth',1.5)
plot(y,t_skin_round(1,:),'b-.','LineWidth',1.5)
plot(y,t_skin_round(10+1,:),'r-.','LineWidth',1.5)
plot(y,t_skin_round(25+1,:),'g-.','LineWidth',1.5)
plot(y,t_skin_round(40+1,:),'k-.','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Skin thickness [m]','Interpreter', 'latex', 'FontSize',18)
legend('0 stringers','10 stringers','25 stringers','40 stringers','Interpreter', 'latex', 'FontSize',16)

for i = 1:length(n_stringers)
    A_stringer(i) = h*(1 + 2*d_h)*ts_t*t_skin(i,1);
    A_stringer_round(i) = h*(1 + 2*d_h)*ts_t*t_skin_round(i,1);
end

As_bt = A_stringer./(b_stringers.*t_skin(:,1));
As_bt_round = A_stringer_round./(b_stringers.*t_skin_round(:,1));

for i = 1:length(y)
    for j = 1:length(n_stringers)
        t_eff(j,i) = t_skin(j,i) + A_stringer(j)./b_stringers(j);
        t_eff_round(j,i) = t_skin_round(j,i) + A_stringer_round(j)./b_stringers(j);
        A_eff(j,i) = t_skin(j,i).*l_box(i) + A_stringer(j).*(floor(l_box(i)./b_stringers(j)));
        A_eff_round(j,i) = t_skin_round(j,i).*l_box(i) + A_stringer_round(j).*(floor(l_box(i)./b_stringers(j)));
    end
end

for i = 1:length(n_stringers)
    V_panel(i) = 0;
    V_panel_round(i) = 0;
    for j = 1:length(y)-1
        V_panel(i) = V_panel(i) + (A_eff(i,j) + A_eff(i,j+1))*0.5*y(2);
        V_panel_round(i) = V_panel_round(i) + (A_eff_round(i,j) + A_eff_round(i,j+1))*0.5*y(2);
    end
end

m_panel = V_panel*rho;
m_panel_round = V_panel_round*rho;

figure
plot(n_stringers,m_panel_CT,'LineWidth',1.5)
hold on
plot(n_stringers,m_panel,'LineWidth',1.5)
plot(n_stringers,m_panel_round,'LineWidth',1.5)
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Number of Stringers','Interpreter', 'latex', 'FontSize',18)
ylabel('Mass [kg]','Interpreter', 'latex', 'FontSize',18)
legend('$t_{skin}$ constant','$t_{skin}$ varying (linearly)','$t_{skin}$ varying (steps)','Interpreter', 'latex', 'FontSize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ribs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rib spacing

for i = 1:length(n_stringers)
    sigma_crit_CT(i) = f1(h_b_CT(i))*E*(t_skin_CT(i)/b_stringers(i)).^2;  % constant thickness
    sigma_crit(i,:) = f0_5(h_b(i))*E*(t_skin(i,:)/b_stringers(i)).^2;
    sigma_crit_round(i,:) = f0_5(h_b(i))*E*(t_skin_round(i,:)/b_stringers(i)).^2;
    sigma_mean_CT(i,:) = (N./t_skin_CT(i))*f1(h_b_CT(i))/3.62;

    Farrar_CT = 0.93;
    Farrar = 0.8;

    L_ribs_CT_final(i,:) = (Farrar_CT./sigma_mean_CT(i,:)).^2.*N*E;
    %L_ribs_CT(i) = min(L_ribs_CT_final(i,:));

    Rb_CT = h_b_CT;
    Rt_CT = ts_t_CT;
    rgyration_CT(i) = sqrt(((b_stringers(i)^2)*Rt_CT*Rb_CT(i)^3*(4+Rt_CT*Rb_CT(i)))/(12*(1+Rt_CT*Rb_CT(i))^2));

    L_ribs_CT_new(i) = pi*b_stringers(i)*rgyration_CT(i)./(t_skin_CT(i)*sqrt(f1(h_b_CT(i))));
    %L_ribs_CT_const(i) = min(L_ribs_CT_new(i));

    Rb = h_b;
    Rt = ts_t;
    rgyration(i) = sqrt(((b_stringers(i)^2)*Rt*Rb(i)^3*(4+Rt*Rb(i)))/(12*(1+Rt*Rb(i))^2));

    L_ribs_new(i,:) = pi*b_stringers(i)*rgyration(i)./(t_skin(i,:)*sqrt(f0_5(h_b(i))));
    L_ribs_const(i) = min(L_ribs_new(i,:));

    L_ribs_round_new(i,:) = pi*b_stringers(i)*rgyration(i)./(t_skin_round(i,:)*sqrt(f0_5(h_b(i))));
    L_ribs_round_const(i) = min(L_ribs_round_new(i,:));
end

% figure
% plot(y,L_ribs_new(1,:),'LineWidth',1.5)
% hold on
% plot(y,L_ribs_new(11,:),'LineWidth',1.5)
% plot(y,L_ribs_new(21,:),'LineWidth',1.5)
% grid on
% box on
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex')
% xlabel('Spanwise Position [m]','Interpreter', 'latex', 'FontSize',18)
% ylabel('Rib spacing [m]','Interpreter', 'latex', 'FontSize',18)
% legend('0 stringers','10 stringers','20 stringers','Interpreter', 'latex', 'FontSize',16)

% L_ribs = transpose((Farrar./sigma_crit(:,1)).^2*N(1)*E);
% L_ribs_round = transpose((Farrar./sigma_crit_round(:,1)).^2*N(1)*E);

L_ribs_CT = L_ribs_CT_new;
L_ribs = L_ribs_const;
L_ribs_round = L_ribs_round_const;

figure
plot(n_stringers,L_ribs_CT,'LineWidth',1.5)
hold on
plot(n_stringers,L_ribs,'LineWidth',1.5)
plot(n_stringers,L_ribs_round,'LineWidth',1.5)
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Number of stringers','Interpreter', 'latex', 'FontSize',18)
ylabel('Rib spacing [m]','Interpreter', 'latex', 'FontSize',18)
legend('$t_{skin}$ constant','$t_{skin}$ varying (linearly)','$t_{skin}$ varying (steps)','Interpreter', 'latex', 'FontSize',16)


% Rib thickness

for i = 1:length(n_stringers)
    A_eff_new(i,:) = t_eff(i,:).*b_stringers(i);
    A_eff_round_new(i,:) = t_eff_round(i,:).*b_stringers(i);
    A_eff_CT_new(i) = t_eff_CT(i).*b_stringers(i);
end

for i = 1:length(n_stringers)
    I(i,:) = (b_stringers(i).*t_eff(i,:).^3)/12 + A_eff_new(i,:).*y_avg/2;
    F(i,:) = (BM_total).^2.*L_ribs(i).*y_avg.*t_eff(i,:).*c./(2*E*I(i,:).^2);
    t_ribs(i,:) = ((F(i,:).*y_avg.^2)./(3.62.*E.*c)).^(1/3);

    I_round(i,:) = (b_stringers(i).*t_eff_round(i,:).^3)/12 + A_eff_round_new(i,:).*y_avg/2;
    F_round(i,:) = (BM_total).^2.*L_ribs_round(i).*y_avg.*t_eff_round(i,:).*c./(2*E*I_round(i,:).^2);
    t_ribs_round(i,:) = ((F_round(i,:).*y_avg.^2)./(3.62.*E.*c)).^(1/3);

    I_CT(i,:) = (b_stringers(i).*t_eff_CT(i).^3)/12 + A_eff_CT_new(i).*y_avg/2;
    F_CT(i,:) = (BM_total).^2.*L_ribs_CT(i).*y_avg.*t_eff_CT(i).*c./(2*E*I_CT(i,:).^2);
    t_ribs_CT(i,:) = ((F_CT(i,:).*y_avg.^2)./(3.62.*E.*c)).^(1/3);
end

t_ribs(t_ribs < 0.001) = 0.001;
t_ribs = ceil(t_ribs*1000)/1000;

t_ribs_CT(t_ribs_round < 0.001) = 0.001;
t_ribs_round = ceil(t_ribs_round*1000)/1000;

t_ribs_CT(t_ribs_CT < 0.001) = 0.001;
t_ribs_CT = ceil(t_ribs_CT*1000)/1000;

figure
plot(y,t_ribs_round(1,:),'b-','LineWidth',1.5)
hold on
plot(y,t_ribs_round(20+1,:),'r-','LineWidth',1.5)
plot(y,t_ribs_round(40+1,:),'g-','LineWidth',1.5)
plot(y,t_ribs_CT(1,:),'b-.','LineWidth',1.5)
plot(y,t_ribs_CT(20+1,:),'r-.','LineWidth',1.5)
plot(y,t_ribs_CT(40+1,:),'g-.','LineWidth',1.5)
plot(y,t_ribs(1,:),'b:','LineWidth',1.5)
plot(y,t_ribs(20+1,:),'r:','LineWidth',1.5)
plot(y,t_ribs(40+1,:),'g:','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise Position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Rib Thickness [m]','Interpreter', 'latex', 'FontSize',18)
legend('0 stringers (SVST)','20 stringers (SVST)','40 stringers (SVST)','0 stringers (CST)','20 stringers (CST)','40 stringers (CST)','0 stringers (LVST)','20 stringers (LVST)','40 stringers (LVST)','Interpreter', 'latex', 'FontSize',16)

% Max number of ribs

Max_ribs_CT = ceil(s./L_ribs_CT);
Max_ribs = ceil(s./L_ribs);
Max_ribs_round = ceil(s./L_ribs_round);

figure
plot(n_stringers,Max_ribs_CT,'b-','LineWidth',1.5)
%hold on
%plot(n_stringers,Max_ribs,'r-','LineWidth',1.5)
%plot(n_stringers,Max_ribs_round,'r-','LineWidth',1.5)
%hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Number of Stringers','Interpreter', 'latex', 'FontSize',18)
ylabel('Number of Ribs','Interpreter', 'latex', 'FontSize',18)
legend('$t_{skin}$ constant','$t_{skin}$ varying (linearly)','$t_{skin}$ varying (steps)','Interpreter', 'latex', 'FontSize',16)

RibSpacingMultiple = floor((length(y)+1)./Max_ribs_CT);
for i = 1:length(n_stringers)
    NextRib = 0;
    j = 1; 
    while NextRib < length(y)
        NextRib = RibSpacingMultiple(i)*j;
        RibLocations(i,j) = NextRib;
        j = j+1;
    end
    NextRib=0;
end
RibLocations(RibLocations > length(y)) = 0;

V_ribs_indiv_CT = zeros(length(n_stringers),size(RibLocations,2));
for i = 1:length(n_stringers)
    for j = 1:size(RibLocations,2)
        if RibLocations(i,j) ~= 0
            V_ribs_indiv_CT(i,j) = l_box(RibLocations(i,j))*y_avg(RibLocations(i,j))*t_ribs_CT(i,RibLocations(i,j));
        end
    end
    V_ribs_CT(i) = sum(V_ribs_indiv_CT(i,:));
end

m_ribs_CT = V_ribs_CT*rho;

V_total_CT = 2*V_panel_CT + V_ribs_CT;

m_total_CT = V_total_CT*rho;


RibSpacingMultiple2 = floor((length(y) + 1)./Max_ribs);
for i = 1:length(n_stringers)
    NextRib = 0;
    j = 1; 
    while NextRib < length(y)
        NextRib = RibSpacingMultiple2(i)*j;
        RibLocations2(i,j) = NextRib;
        j = j+1;
    end
    NextRib=0;
end
RibLocations2(RibLocations2 > length(y)) = 0;

V_ribs_indiv = zeros(length(n_stringers),size(RibLocations2,2));
for i = 1:length(n_stringers)
    for j = 1:size(RibLocations2,2)
        if RibLocations2(i,j) ~= 0
            V_ribs_indiv(i,j) = l_box(RibLocations2(i,j))*y_avg(RibLocations2(i,j))*t_ribs(i,RibLocations2(i,j));
        end
    end
    V_ribs(i) = sum(V_ribs_indiv(i,:));
end

m_ribs = V_ribs*rho;

V_total = 2*V_panel + V_ribs;

m_total = V_total*rho;


RibSpacingMultiple3 = floor((length(y)+1)./Max_ribs_round);
for i = 1:length(n_stringers)
    NextRib = 0;
    j = 1; 
    while NextRib < length(y)
        NextRib = RibSpacingMultiple3(i)*j;
        RibLocations3(i,j) = NextRib;
        j = j+1;
        if NextRib < length(y)
            t_skin_round(i,1+(j-1)*RibSpacingMultiple3(i):NextRib) =  t_skin_round(i,1+(j-1)*RibSpacingMultiple3(i));
        end
    end
    NextRib=0;
end
RibLocations3(RibLocations3 > length(y)) = 0;

V_ribs_indiv_round = zeros(length(n_stringers),size(RibLocations3,2));
for i = 1:length(n_stringers)
    for j = 1:size(RibLocations3,2)
        if RibLocations3(i,j) ~= 0
            V_ribs_indiv_round(i,j) = l_box(RibLocations3(i,j))*y_avg(RibLocations3(i,j))*t_ribs_round(i,RibLocations3(i,j));
        end
    end
    V_ribs_round(i) = sum(V_ribs_indiv_round(i,:));
end

m_panel_round_new = V_panel_round*rho;
m_ribs_round = V_ribs_round*rho;
V_total_round = 2*V_panel_round + V_ribs_round;
m_total_round = V_total_round*rho;

% figure
% plot(y,t_skin(1,:),'b-','LineWidth',1.5)
% hold on
% plot(y,t_skin(5+1,:),'y-','LineWidth',1.5)
% plot(y,t_skin(10+1,:),'r-','LineWidth',1.5)
% plot(y,t_skin(25+1,:),'g-','LineWidth',1.5)
% plot(y,t_skin(50+1,:),'k-','LineWidth',1.5)
% plot(y,t_skin_round(1,:),'b-.','LineWidth',1.5)
% plot(y,t_skin_round(5+1,:),'y-.','LineWidth',1.5)
% plot(y,t_skin_round(10+1,:),'r-.','LineWidth',1.5)
% plot(y,t_skin_round(25+1,:),'g-.','LineWidth',1.5)
% plot(y,t_skin_round(50+1,:),'k-.','LineWidth',1.5)
% hold off
% grid on
% box on
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex')
% xlabel('spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
% ylabel('Skin thickness [m]','Interpreter', 'latex', 'FontSize',18)
% legend('0 stringers','5 stringers','10 stringers','25 stringers','50 stringers','Interpreter', 'latex', 'FontSize',16)

figure
plot(n_stringers,m_total_CT,'b-','LineWidth',1.5)
hold on
plot(n_stringers,m_panel_CT,'r-','LineWidth',1.5)
plot(n_stringers,m_ribs_CT,'g-','LineWidth',1.5)
% plot(n_stringers,m_total_round,'b-.','LineWidth',1.5)
% plot(n_stringers,m_panel_round,'r-.','LineWidth',1.5)
% plot(n_stringers,m_ribs_round,'g-.','LineWidth',1.5)
plot(n_stringers(26),m_total_CT(26),'ko','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Number of Stringers','Interpreter', 'latex', 'FontSize',18)
ylabel('Mass [kg]','Interpreter', 'latex', 'FontSize',18)
legend('Overall','Panel','Ribs','Minimum Point','Interpreter', 'latex', 'FontSize',16)

%% Spars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spars %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q_shear = (SF_total_land)./(2*y_avg);
q_torque = Torque_total_land./(2.*l_box.*y_avg);
q_front = abs(q_shear + q_torque);
q_rear = abs(q_shear - q_torque);

K_s = 8.12;

% Front spar design

t_front_web = ((q_front.*y_avg.^2)./(K_s*E)).^(1/3);

t_front_flange = 0.001:0.001:0.015;
Ixx_front_required = (BM_total.*y_FS/2)/sigma_yield;

for i = 1:length(t_front_flange)
    b_front_flange(i,:) = (Ixx_front_required - (1/12).*t_front_web.*(y_FS - 2*t_front_flange(i)).^3)./((1/6)*t_front_flange(i)^3 + (1/2)*t_front_flange(i).*(y_FS - t_front_flange(i)).^2);
    b_front_flange(b_front_flange < 0.001) = 0.001;
end

for i = 1:length(t_front_flange)
    V_front_flange(i) = 0;
    for j = 1:length(y)-1
        V_front_flange(i) = V_front_flange(i) + (b_front_flange(i,j) + b_front_flange(i,j+1))*t_front_flange(i)*0.5*y(2);
    end
end
m_front_flange = rho*V_front_flange;


% Rear spar design

t_rear_web = ((q_rear.*y_avg.^2)./(K_s*E)).^(1/3);

t_rear_flange = 0.001:0.001:0.015;
Ixx_rear_required = (BM_total.*y_RS/2)/sigma_yield;

for i = 1:length(t_rear_flange)
    b_rear_flange(i,:) = (Ixx_rear_required - (1/12).*t_rear_web.*(y_RS - 2*t_rear_flange(i)).^3)./((1/6)*t_rear_flange(i)^3 + (1/2)*t_rear_flange(i).*(y_RS - t_rear_flange(i)).^2);
    b_rear_flange(b_rear_flange < 0.001) = 0.001;
end

for i = 1:length(t_rear_flange)
    V_rear_flange(i) = 0;
    for j = 1:length(y)-1
        V_rear_flange(i) = V_rear_flange(i) + (b_rear_flange(i,j) + b_rear_flange(i,j+1))*t_rear_flange(i)*0.5*y(2);
    end
end
m_rear_flange = rho*V_rear_flange;


figure
plot(y,t_front_web,'LineWidth',1.5)
hold on
plot(y,t_rear_web,'LineWidth',1.5)
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise Position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Thickness [m]','Interpreter', 'latex', 'FontSize',18)
legend('Front Web','Rear Web','Interpreter', 'latex', 'FontSize',16)

figure
plot(t_front_flange,m_front_flange,'LineWidth',1.5)
hold on
plot(t_rear_flange,m_rear_flange,'LineWidth',1.5)
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Thickness [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Mass [kg]','Interpreter', 'latex', 'FontSize',18)
legend('Front Flange','Rear Flange','Interpreter', 'latex', 'FontSize',16)


%% D-Section

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D-Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% figure
% plot(x1,y1_1,'-b','LineWidth',1.5);
% hold on
% plot(x1,y1_15,'-r','LineWidth',1.5);
% plot(x1,y1_5,'-g','LineWidth',1.5);
% hold off
% box on
% grid on
% xlim([0 20])
% ylim([0 60])
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex')
% xlabel('$a/sqrt(Rt)$','Interpreter', 'latex', 'FontSize',18)
% ylabel('$K_s$','Interpreter', 'latex', 'FontSize',18)
% legend('$1.0$', '$1.5$', '$>5$','Location','southeast','Interpreter', 'latex', 'FontSize',16);
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
% figure
% plot(x2,y2_1,'-b','LineWidth',1.5);
% hold on
% plot(x2,y2_15,'-r','LineWidth',1.5);
% plot(x2,y2_2,'-g','LineWidth',1.5);
% plot(x2,y2_3,'-m','LineWidth',1.5);
% plot(x2,y2_i,'-c','LineWidth',1.5);
% hold off
% grid on
% box on
% xlim([0 20])
% ylim([0 60])
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex')
% xlabel('$s/sqrt(Rt)$','Interpreter', 'latex', 'FontSize',18)
% ylabel('$K_s$','Interpreter', 'latex', 'FontSize',18)
% legend('1.0', '1.5', '2.0', '3.0', 'inf','Location','southeast','Interpreter', 'latex', 'FontSize',16);
% 
% % Save polynomial parameters
% ratiopoly = [p1_1; p1_15; p1_5; p2_1; p2_15; p2_2; p2_3; p2_i];
% 
% tau_tresca = 125*10^6; %%%%%%%%%%%%%%%%%%%%% Update!
% 
% figure
% 
% for k = 1:5 % iterates the whole dcell calculation 3 times to to get converged results
%     n_pseudo_ribs = 0:15;
%     a_pseudo_ribs = s./((n_pseudo_ribs + 1)*cosd(sweep_LE));
%     t_range = [0.001,0.0015,0.002,0.0025,0.003];
%     t_dcell = t_range(k)*ones(1,length(n_pseudo_ribs));
% 
%     data = [0.14645,  0.03994;0.12408,  0.03795;0.10332,  0.03564;0.08427,  0.03305;0.06699,  0.03023;0.05156,  0.02720;0.03806,  0.02395;0.02653,  0.02039;0.01704,  0.01646;0.00961,  0.01214;0.00428,  0.00767;0.00107,  0.00349;0.0,      0.0];
% 
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
% plot([0 15],[max_shear_flow(1),max_shear_flow(1)],'k--','LineWidth',1.5)
% hold off
% grid on
% box on
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex')
% xlabel('Number of Pseudo Ribs','Interpreter', 'latex', 'FontSize',18)
% ylabel('Minimum Buckling Shear Flow [N/m]','Interpreter', 'latex', 'FontSize',18)
% legend('$t = 1\,mm$','$t = 1.5\,mm$','$t = 2\,mm$','$t = 2.5\,mm$','$t = 3\,mm$','Max Applied Shear Flow','Interpreter', 'latex', 'FontSize',18)
% 
% m_dcell = rho*t_range*a_pseudo_ribs(1)*(l_curved(1)+l_curved(end))/2;
% t_pseudo_rib = 0.001;
% 
% for i = 1:length(n_pseudo_ribs)
%     m_pseudo_ribs(i) = 0;
%     if n_pseudo_ribs(i) > 0
%         for j = 1:n_pseudo_ribs(i)
%             m_pseudo_ribs(i) = m_pseudo_ribs(i) + rho*A_pseudo_rib(floor((length(y)*a_pseudo_ribs(i)*j)/(s/cosd(sweep_LE))))*t_pseudo_rib*n_pseudo_ribs(i);
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
% legend('$t = 1\,mm$','$t = 1.5\,mm$','$t = 2\,mm$','$t = 2.5\,mm$','$t = 3\,mm$','Interpreter', 'latex', 'FontSize',18)


%% Wing Configuration Plot

n = 0; % number of stringers
n_pseudo = 0;

ribspacing = [0.204,0.4085,0.4388,0.4388,0.4739,0.4739,0.5642,0.3044];
ribLoc2 = 0;
for i = 1:length(ribspacing)
    ribLoc2(i+1) = ribspacing(i) + ribLoc2(i);
end
figure
p1 = plot([0,s,s,0,0],[c_root,c_root-s*tand(sweep_LE),c_root-s*tand(sweep_LE)-c_tip,0,c_root],'k-','LineWidth',1.5);
hold on 
plot([0,s],[c_root*0.4,c_tip*0.4+c_root-s*tand(sweep_LE)-c_tip],'b-','LineWidth',1.5)
z_LE = @(x) ((c_root-(c_root-s*tand(sweep_LE)))/(0-s))*x + c_root;
z = @(x) ((c_root*0.85-(c_tip*0.85+c_root-s*tand(sweep_LE)-c_tip))/(0-s))*x + c_root*0.85;
z_RS = @(x) ((c_root*0.4-(c_tip*0.4+c_root-s*tand(sweep_LE)-c_tip))/(0-s))*x + c_root*0.4;
for i = 1:n
    H = c_tip*0.4+c_root-s*tand(sweep_LE)-c_tip+i*b_stringers(n+1);
    if H <= c_tip*0.85+c_root-s*tand(sweep_LE)-c_tip
        p3 = plot([0,s],[c_root*0.4+i*b_stringers(n+1),H],'g-','LineWidth',1.5);
    else
        z_stringer = @(x) ((c_root*0.4+i*b_stringers(n+1)-H)/(0-s))*x + c_root*0.4+i*b_stringers(n+1);
        Int = fzero(@(x) z(x)-z_stringer(x), 1);
        plot([0,Int],[c_root*0.4+i*b_stringers(n+1),z(Int)],'g-','LineWidth',1.5)
    end
end
for j = 2:length(ribLoc2)-1
    p4 = plot([0+ribLoc2(j),0+ribLoc2(j)],[z(ribLoc2(j)),z_RS(ribLoc2(j))],'r-','LineWidth',1.5); %z(ribLoc2(j))-l_box(round(ribLoc2(j)*100/s)
end
if n_pseudo > 0
    for k = 1:n_pseudo
        spacing = floor(length(y)/(n_pseudo+1));
        p5 = plot([0+y(spacing*k),0+y(spacing*k)],[z_LE(y(spacing*k)),z(y(spacing*k))],'y-','LineWidth',1.5);
    end
end
p2 = plot([0,s],[c_root*0.85,c_tip*0.85+c_root-s*tand(sweep_LE)-c_tip],'b-','LineWidth',1.5);
p6 = fill([0,s*0.8,s*0.8,0,0],[c_root*0.25,c_root*0.25-s*0.8*tand(sweep_hinge),c_root*0.25-s*0.8*tand(sweep_hinge)-0.3,0,c_root*0.25],[0.827,0.827,0.827]); %0 0.9 0.5
hold off
box on 
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise Position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Chordwise Position [m]','Interpreter', 'latex', 'FontSize',18)
legend([p2 p4 p6],{'Spars','Ribs','Elevators'},'Interpreter', 'latex', 'FontSize',16)

