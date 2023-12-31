clear
clc

%Load Case 1

n = 3.75;
n_low = -1.5;
n_A = 2.5;
n_A_low = -1;

y_LC1 = [0, 0.25, 0.5, 1, 1.5, 1.67, 2, 2.5, 3, 3.5, 3.74, 3.94, 4, 4.14, 4.5, 5, 5.5, 5.86, 6, 6.5, 6.56, 7, 7.25, 7.5, 7.95, 8, 8.5, 8.65, 9, 9.35, 9.37, 9.5, 10, 10.05, 10.5, 10.75, 11, 11.44, 11.5, 11.93, 12, 12.14, 12.5, 12.84, 13, 13.3, 13.5, 13.71, 14, 14.5, 14.75, 15, 15.11, 15.46, 15.5, 16, 16.16, 16.5, 16.86, 16.93, 17, 17.5, 17.57, 17.75, 18, 18.27, 18.5, 18.62, 18.97, 19, 19.5, 19.68, 20, 20.38, 20.5, 21, 21.08, 21.5, 21.79, 22, 22.5, 23, 23.19, 23.5, 24, 24.41, 24.5, 24.61, 25, 25.25, 25.5, 26, 26.5, 26.81, 27, 27.5, 28, 28.5, 29, 29.49, 29.5, 29.656, 30, 30.5, 31, 31.5];

W = [-612.3635571, -1017.5913, -612.3635571, -6109.298957, -612.3635571, -2358.6183, -2961.368057, -612.3635571, -612.3635571, -612.3635571, -4905, -1948.8546, -891.9485571, -3653.6364, -891.9485571, -891.9485571, -891.9485571, -2109.15, -891.9485571, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -891.9485571, -4218.3, -891.9485571, -4218.3, -4905, -891.9485571, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -1002.4839, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, 28300.94159, -891.9485571, 0, -891.9485571, -891.9485571, -4218.3, -612.3635571, 192147.1635, -4218.3, -612.3635571, -612.3635571, -4218.3, -891.9485571, -4218.3, -3626.0703, -891.9485571, -891.9485571, -4218.3, -12303.0153, -891.9485571, -4218.3, -891.9485571, -4218.3, -4218.3, -891.9485571, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -891.9485571, -891.9485571, -12148.5078, -891.9485571, -891.9485571, -1948.8546, -612.3635571, -4905, -612.3635571, -10837.9899, -612.3635571, -612.3635571, -612.3635571, -707.1048, -612.3635571, -612.3635571, -612.3635571, -612.3635571, -612.3635571, -1193.1903, -612.3635571, -600, -612.3635571, -612.3635571, -612.3635571, -612.3635571];

W_LC1 = n*W;
W_LC1_low = n_low*W;
W_LC1_A = n_A*W;
W_LC1_A_low = n_A_low*W;

%Load Case 3 (Flare Landing)


%Calculation
% y = [0, 0.25, 0.5, 1, 1.67, 1.5, 2, 2.5, 3, 3.5, 4.14, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 23.73, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28, 28.5, 29, 29.5, 30, 30.5, 31, 31.5];
% 
% M_inertia = [62.42, 103.73, 62.42, 560.34, 62.42, 62.42, 186.00, 54.43, 239.45, 62.42, 62.42, 62.42, 62.42, 500.00, 198.66, 90.92, 279.00, 81.65, 11.79, 90.92, 90.92, 90.92, 215.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 430.00, 90.92, 430.00, 0.00, 500.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 430.00, 90.92, 102.19, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 0.00, 90.92, 90.92, 430.00, 120.75, 7601.00, 62.42, 2856.00, 4402.17, 125.96, 78.58, 380.00, 430.00, 62.42, 62.42, 430.00, 0.00, 90.92, 430.00, 369.63, 516.70, 90.92, 90.92, 430.00, 1254.13, 90.92, 430.00, 90.92, 430.00, 0.00, 430.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 90.92, 1238.38, 90.92, 90.92, 198.66, 62.42, 500.00, 62.42, 1104.79, 62.42, 62.42, 62.42, 72.08, 62.42, 62.42, 62.42, 62.42, 62.42, 121.63, 62.42, 62.42, 62.42, 62.42, 62.42];
% 
% M_total = sum(M_inertia);
% 
% MTOW = M_total*9.81; % [N]
% n = 3.75;
% Sref_h = 10;
% taper_h = 0.45;
% AR = 4;
% sweep = 17;
% b = sqrt(AR*Sref_h);
% c_root = 2*Sref_h/(b*(1+taper_h));
% 
% 
% x_h = 28.739 + 0.25*c_root;
% x_cg = 15.08;
% x_mlg = 17.71;
% 
% 
% % F_h + F_mlg = MTOW (sum of vertical forces = 0)
% 
% % F_h*(x_h - x_cg) + F_mlg*(x_mlg - x_cg) = 0 (sum of moments = 0)
% 
% % Solving simultaneous equation gives:
% 
% F_h = -MTOW*(x_mlg - x_cg)/((x_h - x_cg) - (x_mlg - x_cg));
% F_mlg = MTOW - F_h;
% 
% y_land_flare = [0, 0.25, 0.5, 1, 1.5, 1.67, 2, 2.5, 3, 3.5, 3.74, 3.94, 4, 4.14, 4.5, 5, 5.5, 5.86, 6, 6.5, 6.56, 7, 7.25, 7.5, 7.95, 8, 8.5, 8.65, 9, 9.35, 9.37, 9.5, 10, 10.05, 10.5, 10.75, 11, 11.44, 11.5, 11.93, 12, 12.14, 12.5, 12.84, 13, 13.5, 13.71, 14, 14.5, 14.75, 14.86, 14.86, 15, 15.02, 15.14, 15.14, 15.14, 15.15, 15.46, 15.5, 16, 16.16, 16.5, 16.86, 16.93, 16.96, 17, 17.5, 17.57, 17.71, 17.75, 18, 18.27, 18.5, 18.62, 18.97, 19, 19.5, 19.68, 20, 20.38, 20.5, 21, 21.08, 21.5, 21.79, 22, 22.5, 23, 23.19, 23.5, 24, 24.41, 24.5, 24.61, 25, 25.25, 25.5, 26, 26.5, 26.81, 27, 27.5, 28, 28.5, 29, 29.284, 29.49, 29.5, 30, 30.5, 31, 31.5];
% 
% W_inertia_land_flare = 9.81*[62.42, 103.73, 62.42, 560.34+62.42, 62.42, 186.00+54.43, 239.45+62.42, 62.42, 62.42, 62.42, 500.00, 198.66, 90.92, 279.00+81.65+11.79, 90.92, 90.92, 90.92, 215.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 430.00, 90.92, 430.00, 500.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 430.00, 90.92, 102.19, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 90.92, 90.92, 430.00, 120.75, 7475.00, 62.42, 2856.00, 4402.17, 125.96, 78.58, 380.00, 430.00, 62.42, 62.42, 430.00, 90.92, 430.00, 369.63, 516.70, 90.92, 90.92, 430.00, -F_mlg/9.81, 1254.13, 90.92, 430.00, 90.92, 430.00, 0.00, 430.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 90.92, 1238.38, 90.92, 90.92, 198.66, 62.42, 500.00, 62.42, 1104.79, 62.42, 62.42, 62.42, 72.08, 62.42, 62.42, 62.42, 62.42, 62.42, -F_h/9.81, 121.63, 62.42, 62.42, 62.42, 62.42, 62.42];

y_land = [0, 0.25, 0.5, 1, 1.5, 1.67, 2, 2.5, 3, 3.5, 3.74, 3.94, 4, 4.14, 4.5, 5, 5.5, 5.86, 6, 6.5, 6.56, 7, 7.25, 7.5, 7.95, 8, 8.5, 8.65, 9, 9.35, 9.37, 9.5, 10, 10.05, 10.5, 10.75, 11, 11.44, 11.5, 11.93, 12, 12.14, 12.5, 12.84, 13, 13.3, 13.5, 13.71, 14, 14.5, 14.75, 15, 15.11, 15.46, 15.5, 16, 16.16, 16.5, 16.86, 16.93, 17, 17.5, 17.57, 17.71, 17.75, 18, 18.27, 18.5, 18.62, 18.97, 19, 19.5, 19.68, 20, 20.38, 20.5, 21, 21.08, 21.5, 21.79, 22, 22.5, 23, 23.19, 23.5, 24, 24.41, 24.5, 24.61, 25, 25.25, 25.5, 26, 26.5, 26.81, 27, 27.5, 28, 28.5, 29, 29.49, 29.5, 29.656, 30, 30.5, 31, 31.5];

W_land = [-612.3635571, -1017.5913, -612.3635571, -6109.298957, -612.3635571, -2358.6183, -2961.368057, -612.3635571, -612.3635571, -612.3635571, -4905, -1948.8546, -891.9485571, -3653.6364, -891.9485571, -891.9485571, -891.9485571, -2109.15, -891.9485571, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -891.9485571, -4218.3, -891.9485571, -4218.3, -4905, -891.9485571, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -1002.4839, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -16766.87044, -891.9485571, 0, -891.9485571, -891.9485571, -4218.3, -612.3635571, -424.511963, -4218.3, -612.3635571, -612.3635571, -4218.3, -891.9485571, -4218.3, -3626.0703, -891.9485571, -891.9485571, -4218.3, 295639, -12303.0153, -891.9485571, -4218.3, -891.9485571, -4218.3, -4218.3, -891.9485571, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -891.9485571, -4218.3, -891.9485571, -4218.3, -891.9485571, -891.9485571, -891.9485571, -12148.5078, -891.9485571, -891.9485571, -1948.8546, -612.3635571, -4905, -612.3635571, -10837.9899, -612.3635571, -612.3635571, -612.3635571, -707.1048, -612.3635571, -612.3635571, -612.3635571, -612.3635571, -612.3635571, -1193.1903, -612.3635571, -58413.5, -612.3635571, -612.3635571, -612.3635571, -612.3635571];


% Load Case 3 (Level Landing)

% x_nlg = 2;
% x_m = x_mlg - x_nlg;
% x_g = x_cg - x_nlg;
% 
% F_m = -x_g*(-MTOW)/x_m;
% F_n = MTOW-F_m;
% disp(F_m)
% disp(F_n)
% 
% F_mlg2 = -MTOW*(x_nlg - x_cg)/((x_mlg - x_cg) - (x_nlg - x_cg));
% F_nlg = MTOW - F_mlg2;
% 
% y_land_level = [0, 0.25, 0.5, 1, 1, 1.5, 1.67, 1.67, 2, 2, 2, 2.5, 3, 3.5, 3.74, 3.94, 4, 4.14, 4.14, 4.14, 4.5, 5, 5.5, 5.86, 6, 6.5, 6.56, 7, 7.25, 7.5, 7.95, 8, 8.5, 8.65, 9, 9.35, 9.37, 9.5, 10, 10.05, 10.5, 10.75, 11, 11.44, 11.5, 11.93, 12, 12.14, 12.5, 12.84, 13, 13.5, 13.71, 14, 14.5, 14.75, 14.86, 14.86, 15, 15.02, 15.14, 15.14, 15.14, 15.15, 15.46, 15.5, 16, 16.16, 16.5, 16.86, 16.93, 16.96, 17, 17.5, 17.57, 17.71, 17.75, 18, 18.27, 18.5, 18.62, 18.97, 19, 19.5, 19.68, 20, 20.38, 20.5, 21, 21.08, 21.5, 21.79, 22, 22.5, 23, 23.19, 23.5, 24, 24.41, 24.5, 24.61, 25, 25.25, 25.5, 26, 26.5, 26.81, 27, 27.5, 28, 28.5, 29, 29.49, 29.5, 30, 30.5, 31, 31.5];
% 
% W_inertia_land_level = 9.81*[62.42, 103.73, 62.42, 560.34, 62.42, 62.42, 186.00, 54.43, -7177.78, 239.45, 62.42, 62.42, 62.42, 62.42, 500.00, 198.66, 90.92, 279.00, 81.65, 11.79, 90.92, 90.92, 90.92, 215.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 430.00, 90.92, 430.00, 500.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 430.00, 90.92, 430.00, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 3932.61, 90.92, 90.92, 430.00, 120.75, 7475.00, 90.92, 2856.00, 4402.17, 125.96, 78.58, 380.00, 430.00, 62.42, 90.92, 430.00, 90.92, 430.00, 369.63, 516.70, 90.92, 90.92, 430.00, -32112.13, 1254.13, 90.92, 430.00, 90.92, 430.00, 430.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 430.00, 90.92, 430.00, 90.92, 90.92, 90.92, 1238.38, 62.42, 90.92, 198.66, 62.42, 500.00, 62.42, 1104.79, 62.42, 62.42, 62.42, 72.08, 62.42, 62.42, 62.42, 62.42, 62.42, 121.63, 62.42, 62.42, 62.42, 62.42, 62.42];
% 

SF = zeros(1,length(y_LC1));
dM = zeros(1,length(y_LC1));
BM = zeros(1,length(y_LC1));
SF_low = zeros(1,length(y_LC1));
dM_low = zeros(1,length(y_LC1));
BM_low = zeros(1,length(y_LC1));
SF_land = zeros(1,length(y_land));
dM_land = zeros(1,length(y_land));
BM_land = zeros(1,length(y_land));
SF_A = zeros(1,length(y_LC1));
dM_A = zeros(1,length(y_LC1));
BM_A = zeros(1,length(y_LC1));
SF_A_low = zeros(1,length(y_LC1));
dM_A_low = zeros(1,length(y_LC1));
BM_A_low = zeros(1,length(y_LC1));

for i = 1:length(y_LC1)
    SF(i) = sum(W_LC1(1:i));
    SF_low(i) = sum(W_LC1_low(1:i));
    SF_A(i) = sum(W_LC1_A(1:i));
    SF_A_low(i) = sum(W_LC1_A_low(1:i));
end

for i = 2:length(y_LC1)
    dM(i) = (SF(i) + SF(i-1))*(y_LC1(i) - y_LC1(i-1))*0.5;
    dM_low(i) = (SF_low(i) + SF_low(i-1))*(y_LC1(i) - y_LC1(i-1))*0.5;
    dM_A(i) = (SF_A(i) + SF_A(i-1))*(y_LC1(i) - y_LC1(i-1))*0.5;
    dM_A_low(i) = (SF_A_low(i) + SF_A_low(i-1))*(y_LC1(i) - y_LC1(i-1))*0.5;
end

for i = 1:length(y_LC1)
    BM(i) = sum(dM(1:i));
    BM_low(i) = sum(dM_low(1:i));
    BM_A(i) = sum(dM_A(1:i));
    BM_A_low(i) = sum(dM_A_low(1:i));
end


for i = 1:length(y_land)
    SF_land(i) = sum(W_land(1:i));
end

for i = 2:length(y_land)
    dM_land(i) = (SF_land(i) + SF_land(i-1))*(y_land(i) - y_land(i-1))*0.5;
end

for i = 1:length(y_land)
    BM_land(i) = sum(dM_land(1:i));
end


figure
plot(y_LC1,SF,'b-','LineWidth',1.5)
hold on
plot(y_LC1,SF_low,'b-.','LineWidth',1.5)
plot(y_LC1,SF_A,'k-','LineWidth',1.5)
plot(y_LC1,SF_A_low,'k-.','LineWidth',1.5)
plot(y_land,SF_land,'r-','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Shear Force [N]','Interpreter', 'latex', 'FontSize',18)
legend('LC1 ($V_D$ at n = 3.75)','LC1 ($V_D$ at n = -1.5)','LC1 ($V_A$ at n = 2.5)','LC1 ($V_A$ at n = -1)','LC3 (landing)','Interpreter', 'latex', 'FontSize',16)

figure
plot(y_LC1,BM,'b-','LineWidth',1.5)
hold on
plot(y_LC1,BM_low,'b-.','LineWidth',1.5)
plot(y_LC1,BM_A,'k-','LineWidth',1.5)
plot(y_LC1,BM_A_low,'k-.','LineWidth',1.5)
plot(y_land,BM_land,'r-','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Bending Moment [Nm]','Interpreter', 'latex', 'FontSize',18)
legend('LC1 ($V_D$ at n = 3.75)','LC1 ($V_D$ at n = -1.5)','LC1 ($V_A$ at n = 2.5)','LC1 ($V_A$ at n = -1)','LC3 (landing)','Interpreter', 'latex', 'FontSize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Torque %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sref_VT = 12.3;
taper_VT = 0.7;
AR_VT = 0.9;
sweep_VT = 40;
b_VT = sqrt(AR_VT*Sref_VT);
c_root_VT = 2*Sref_VT/(b_VT*(1+taper_VT));
c_tip_VT = taper_VT*c_root_VT; 
c_root_w = 4.0259;

y_bar_VT = (b_VT/6)*((1+2*taper_VT)*(1+taper_VT));

x_VT = 27.17 + 0.25*c_root_VT;
x_w = 12.7;

% Find lift at vertical tailplane for moment balance with OEI

L_VT = 63000*14.1*0.3/(x_VT - (x_w +0.5*c_root_w));

y_VT = 1.343;  %%%% Check as values from engineering drawing and designparams slightly different

T = y_VT*L_VT;

figure
plot([0,x_w,x_w,x_VT,x_VT,31.5],[0,0,T,T,0,0],'g-','LineWidth',1.5)
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Torque [Nm]','Interpreter', 'latex', 'FontSize',18)
legend('Load Case 2','Interpreter', 'latex', 'FontSize',16)


figure
h1 = bar(y_LC1,W_LC1,'b','BarWidth',15);
hold on
h2 = bar(y_land,W_land,'r','BarWidth',15);
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Distributed Load [N/m]','Interpreter', 'latex', 'FontSize',18)
legend('Load Case 1 (n=1)','Load Case 3 (landing)','Interpreter', 'latex', 'FontSize',16)
h1.EdgeColor = 'none';
h2.EdgeColor = 'none';
