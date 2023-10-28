clear
clc

x_FA = 0.35; % Flexural Axis 
x_AC = 0.25; % Aerodynamic Centre
x_CG = 0.4; % Wing Centre of Gravity (update to more accurate one!!)

Sref = 12.3;
taper = 0.7;
AR = 0.9;
sweep_LE = 40;
rho = 1.225;
sweep = atand(tand(sweep_LE)-(4/AR/2)*((x_AC - 0)*(1-taper)/(1+taper))); % at quarter chord
sweep_FA = atand(tand(sweep)-(4/AR/2)*((x_FA - x_AC)*(1-taper)/(1+taper)));
sweep_hinge = atand(tand(sweep_LE)-(4/AR/2)*((0.68 - 0)*(1-taper)/(1+taper)));

b = sqrt(AR*Sref);
s = b;
y = linspace(0,s,1000);
c_root = 2*Sref/(b*(1+taper));
c_tip = taper*c_root; 
c_root_w = 4.0259;

x_VT = 27.17 + 0.25*c_root;
x_w = 12.7;

% Find lift at vertical tailplane for moment balance with OEI

L_VT = 1.5*(63000*1.15*14.1*0.3/(x_VT - (x_w +0.5*c_root_w))); % Safety factor of 1.5 (as there are no load factors for VTP loads)

% In this case lift comes only from one wing

L0_VT = 4*L_VT/(pi*s);
dL_VT = L0_VT*sqrt(1-(y/s).^2);

Lift = zeros(1,length(y));
SF_aero = zeros(1,length(y));
dM_aero = zeros(1,length(y));
BM_aero = zeros(1,length(y));

for i = 1:length(y)-1
    Lift(i) = (dL_VT(i+1) + dL_VT(i))*(y(i+1) - y(i))*0.5;
end

for i = 1:length(y)-1
    SF_aero(i) = sum(Lift(i:end));
end

for i = 1:length(y)-1
    dM_aero(i) = (SF_aero(i+1) + SF_aero(i))*(y(i+1) - y(i))*0.5;
end

for i = 1:length(y)-1
    BM_aero(i) = sum(dM_aero(i:end));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Torque %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = (c_tip - c_root)*(y/s) + c_root;

dM_lift = -dL_VT.*c*(x_FA - x_AC);

M_lift = zeros(1,length(y));

for i = 1:length(y)-1
    M_lift(i) = (dM_lift(i+1) + dM_lift(i))*(y(i+1) - y(i))*0.5;
end

dT = M_lift;

Torque = zeros(1,length(y));

for i = 1:length(y)-1
    Torque(i) = sum(dT(i:end));
end

Torque_total = Torque + BM_aero*sind(sweep_FA);
BM_total = BM_aero + Torque*sind(sweep_FA);


figure
plot(y,dL_VT,'g-','LineWidth',1.5)
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Load Distribution [N/m]','Interpreter', 'latex', 'FontSize',18)
legend('Lift','Interpreter', 'latex', 'FontSize',16)

figure
plot(y,SF_aero,'g-','LineWidth',1.5)
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Shear Force [N]','Interpreter', 'latex', 'FontSize',18)
legend('Load Case 2','Interpreter', 'latex', 'FontSize',16)

figure
plot(y,BM_total,'g-','LineWidth',1.5)
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Bending Moment [Nm]','Interpreter', 'latex', 'FontSize',18)
legend('Load Case 2','Interpreter', 'latex', 'FontSize',16)

figure
plot(y,Torque_total,'g-','LineWidth',1.5)
box on 
grid on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Torque [Nm]','Interpreter', 'latex', 'FontSize',18)
legend('Load Case 2','Interpreter', 'latex', 'FontSize',16)


%% Design 

% Parameters
% Al 7075-T6 Material Properties   Al 2024 T861
E = 73.8*10^9; 
rho = 2.77*100^3/1000;
sigma_yield = 430*10^6;

x_FS = 0.15;
x_RS = 0.55;
y_FS_c1 = 0.08038; % From airfoil tools (For c = 1)
y_RS_c1 = 0.07436; % From airfoil tools (For c = 1)
y_avg_c1 = (y_FS_c1 + y_RS_c1)/2;

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


N_root = abs((BM_total(1))./(l_box(1).*y_avg(1))); %compressive load per unit length


%%%%%%%%%%%%%%%%%%%%% stringers (t_skin constant) %%%%%%%%%%%%%%%%%%%%%%%%%

n_stringers = 0:30;
b_stringers = l_box(1)./(n_stringers + 1);

d_h_CT = 0.3; % Value indicated on ESDU plot

ts_t_CT = 1; % selected for max farrar  efficiency

h = 0.125*y_avg(end);  % 35% of rear spar height at tip
h_b_CT = h./b_stringers;

hb_1 = [0.101663586, 0.12754159, 0.153419593, 0.22181146, 0.273567468, 0.336414048, 0.410351201, 0.46025878, 0.547134935, 0.600739372, 0.641404806, 0.689463956, 0.748613678, 0.802218115, 0.859519409, 0.896487985, 0.926062847, 0.959334566, 1];
K_1 = [5.488946459, 5.383592401, 5.278238342, 5.056994819, 4.920034542, 4.814680484, 4.709326425, 4.677720207, 4.603972366, 4.561830743, 4.519689119, 4.46701209, 4.393264249, 4.308981002, 4.214162349, 4.129879102, 4.066666667, 3.961312608, 3.813816926];

[xData4, yData4] = prepareCurveData(hb_1, K_1);
f1 = fit(xData4, yData4,'pchip');

t_skin_CT = zeros(1,length(n_stringers));


for i = 1:length(n_stringers)
    t_skin_CT(i) = ((N_root.*b_stringers(i).^2)./(3.62*E)).^(1/3);
end
t_skin_CT(t_skin_CT < 0.001) = 0.001;

A_stringer_CT = h*(1 + 2*d_h_CT)*ts_t_CT*t_skin_CT;
As_bt_CT = A_stringer_CT./(b_stringers.*t_skin_CT);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ribs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rib spacing

N = abs((BM_total)./(l_box.*y_avg)); %compressive load per unit length

for i = 1:length(n_stringers)
    sigma_crit_CT(i) = f1(h_b_CT(i))*E*(t_skin_CT(i)/b_stringers(i)).^2;  % constant thickness
    sigma_applied(i,:) = N/t_eff_CT(i);
    sigma_applied_max(i) = max(sigma_applied(i,:));
    
    sqrt_I_A(i) = sqrt(h^2*h_b_CT(i)*ts_t_CT*(0.633+0.37*h_b_CT(i)*ts_t_CT)/((1.6*h_b_CT(i)*ts_t_CT + 1)^2));
    L_ribs_CT_new(i) = pi*b_stringers(i)*sqrt_I_A(i)./(t_skin_CT(i)*sqrt(f1(h_b_CT(i))));
    %L_ribs_CT_const(i) = min(L_ribs_CT_new(i));

end

figure
plot(n_stringers,sigma_applied_max,'b-','LineWidth',1.5)
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Number of Stringers','Interpreter', 'latex', 'FontSize',18)
ylabel('Max Stress [Pa]','Interpreter', 'latex', 'FontSize',18)



L_ribs_CT = L_ribs_CT_new;

% Rib thickness

for i = 1:length(n_stringers)
    I_CT(i,:) = (b_stringers(i).*t_eff_CT(i).^3)/12 + A_eff_CT(i,:).*y_avg/2;
    F_CT(i,:) = (BM_total).^2.*L_ribs_CT(i).*y_avg.*t_eff_CT(i).*c./(2*E*I_CT(i,:).^2);
    t_ribs_CT(i,:) = ((F_CT(i,:).*y_avg.^2)./(3.62.*E.*c)).^(1/3);
end

t_ribs_CT(t_ribs_CT < 0.001) = 0.001;
t_ribs_CT = ceil(t_ribs_CT*1000)/1000;

figure
plot(y,t_ribs_CT(1,:),'b-.','LineWidth',1.5)
hold on 
plot(y,t_ribs_CT(20+1,:),'r-.','LineWidth',1.5)
plot(y,t_ribs_CT(30+1,:),'g-.','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise Position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Rib Thickness [m]','Interpreter', 'latex', 'FontSize',18)
legend('0 stringers (CST)','20 stringers (CST)','30 stringers (CST)','Interpreter', 'latex', 'FontSize',16)



% Max number of ribs

Max_ribs_CT = ceil(s./L_ribs_CT);

figure
plot(n_stringers,Max_ribs_CT,'b-','LineWidth',1.5)
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Number of Stringers','Interpreter', 'latex', 'FontSize',18)
ylabel('Number of Ribs','Interpreter', 'latex', 'FontSize',18)
legend('$t_{skin}$ constant','$t_{skin}$ varying (linearly)','$t_{skin}$ varying (steps)','Interpreter', 'latex', 'FontSize',16)

RibSpacingMultiple = ceil((length(y))./Max_ribs_CT);
for i = 1:length(n_stringers)
    NextRib = 0;
    j = 1;
    if Max_ribs_CT(i) > 0
        while NextRib < length(y)
            NextRib = RibSpacingMultiple(i)*j;
            RibLocations(i,j) = NextRib;
            j = j+1;
        end
    else
        RibLocations(i,j) = 0;
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


figure
plot(n_stringers,m_total_CT,'b-','LineWidth',1.5)
hold on
plot(n_stringers,m_panel_CT,'r-','LineWidth',1.5)
plot(n_stringers,m_ribs_CT,'g-','LineWidth',1.5)
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

h_percent = 0.05:0.025:0.35;
for i = 1:length(h_percent)
    [m_total(i,:),n_ribs(i,:)] = SkinAndRibsMass(h_percent(i));
end

figure
surf(n_stringers,h_percent*100,real(m_total),n_ribs-1)
C = colorbar;
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Number of Stringers','Interpreter', 'latex', 'FontSize',18)
ylabel('Stringer Height [\%]','Interpreter', 'latex', 'FontSize',18)
zlabel('Mass [kg]','Interpreter', 'latex', 'FontSize',18)
ylabel(C,'Number of Ribs','Interpreter', 'latex', 'FontSize',18)
view(40,40);


%% Spars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spars %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 25; % number of stringers

q_shear = (SF_aero)./(2*y_avg);
q_torque = Torque_total./(2.*l_box.*y_avg);
q_front = abs(q_shear + q_torque);
q_rear = abs(q_shear - q_torque);

K_s = 8.12;

% Front spar design

t_front_web = ((q_front.*y_avg.^2)./(K_s*E)).^(1/3);
t_front_web_steps = ceil(t_front_web*10000)/10000;
t_front_web_steps(t_front_web_steps < 0.001) = 0.001;
t_front_web_steps_final = t_front_web_steps;
for j = 1:size(RibLocations,2)
    if RibLocations(i,j) ~= 0 && j == 1
        t_front_web_steps_final(1+(j-1)*RibSpacingMultiple(n+1):end) =  t_front_web_steps(1+(j-1)*RibSpacingMultiple(n+1));
    end
    if RibLocations(i,j) ~= 0 && j > 1
        t_front_web_steps_final((j-1)*RibSpacingMultiple(n+1):end) =  t_front_web_steps((j-1)*RibSpacingMultiple(n+1));
    end
end
%t_front_web_steps((j-1)*RibSpacingMultiple(n+1):end) = t_front_web_steps(1+(j-2)*RibSpacingMultiple(n+1));


t_front_flange = 0.001:0.001:0.015;
Ixx_front_required = (BM_total.*y_avg/2)/sigma_yield;

for i = 1:length(t_front_flange)
    b_front_flange(i,:) = (Ixx_front_required - (1/12).*t_front_web_steps_final.*(y_avg - 2*t_front_flange(i)).^3)./((1/6)*t_front_flange(i)^3 + (1/2)*t_front_flange(i).*(y_avg - t_front_flange(i)).^2);
    b_front_flange_max(i) = max(b_front_flange(i,:));
end
b_front_flange(b_front_flange < 0.01) = 0.01; % 1cm lower bound to allow for riveting

for i = 1:length(t_front_flange)
    V_front_flange(i) = 0;
    V_front_web(i) = 0;
    for j = 1:length(y)-1
        V_front_flange(i) = V_front_flange(i) + (b_front_flange(i,j) + b_front_flange(i,j+1))*t_front_flange(i)*0.5*y(2);
        V_front_web(i) = V_front_web(i) + (y_avg(j) + y_avg(j+1))*t_front_web_steps_final(j)*0.5*y(2);
    end
end
m_front_flange = rho*V_front_flange;
m_front_web = rho*V_front_web;
m_front_spar = m_front_web(1) + m_front_flange(3);


% Rear spar design

t_rear_web = ((q_rear.*y_avg.^2)./(K_s*E)).^(1/3);
t_rear_web_steps = ceil(t_rear_web*10000)/10000;
t_rear_web_steps(t_rear_web_steps < 0.001) = 0.001;
t_rear_web_steps_final = t_rear_web_steps;
for j = 1:size(RibLocations,2)
    if RibLocations(i,j) ~= 0 && j == 1
        t_rear_web_steps_final(1+(j-1)*RibSpacingMultiple(n+1):end) =  t_rear_web_steps(1+(j-1)*RibSpacingMultiple(n+1));
    end
    if RibLocations(i,j) ~= 0 && j > 1
        t_rear_web_steps_final((j-1)*RibSpacingMultiple(n+1):end) =  t_rear_web_steps((j-1)*RibSpacingMultiple(n+1));
    end
end
%t_rear_web_steps((j-1)*RibSpacingMultiple(n+1):end) = t_rear_web_steps(1+(j-2)*RibSpacingMultiple(n+1));

t_rear_flange = 0.001:0.001:0.015;
Ixx_rear_required = (BM_total.*y_avg/2)/sigma_yield;

for i = 1:length(t_rear_flange)
    b_rear_flange(i,:) = (Ixx_rear_required - (1/12).*t_rear_web_steps_final.*(y_avg - 2*t_rear_flange(i)).^3)./((1/6)*t_rear_flange(i)^3 + (1/2)*t_rear_flange(i).*(y_avg - t_rear_flange(i)).^2);
    b_rear_flange_max(i) = max(b_rear_flange(i,:));
end
b_rear_flange(b_rear_flange < 0.01) = 0.01; % 1cm lower bound to allow for riveting

for i = 1:length(t_rear_flange)
    V_rear_flange(i) = 0;
    V_rear_web(i) = 0;
    for j = 1:length(y)-1
        V_rear_flange(i) = V_rear_flange(i) + (b_rear_flange(i,j) + b_rear_flange(i,j+1))*t_rear_flange(i)*0.5*y(2);
        V_rear_web(i) = V_rear_web(i) + (y_avg(j) + y_avg(j+1))*t_rear_web_steps_final(j)*0.5*y(2);
    end
end
m_rear_flange = rho*V_rear_flange;
m_rear_web = rho*V_rear_web;
m_rear_spar = m_rear_flange(3) + m_rear_web(1);

figure
plot(y,t_front_web,'b-','LineWidth',1.5)
hold on
plot(y,t_front_web_steps_final,'b--','LineWidth',1.5)
plot(y,t_rear_web,'r-','LineWidth',1.5)
plot(y,t_rear_web_steps_final,'r--','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise Position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Web Thickness [m]','Interpreter', 'latex', 'FontSize',18)
legend('Front Spar (Ideal)','Front Spar (Actual)','Rear Spar (Ideal)','Rear Spar (Actual)','Interpreter', 'latex', 'FontSize',16)

figure
plot(t_front_flange,b_front_flange_max,'LineWidth',1.5)
hold on
plot(t_rear_flange,b_rear_flange_max,'LineWidth',1.5)
yline(b_stringers(26),'k--','LineWidth',1.5)
plot(0.003,b_front_flange_max(3),'o','Color',[0 0.4470 0.7410],'LineWidth',1.5)
plot([0.003,0.003],[b_front_flange_max(3),0],'-.','Color',[0 0.4470 0.7410],'LineWidth',1.5)
plot(0.003,b_rear_flange_max(3),'o','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
plot([0.003,0.003],[b_rear_flange_max(3),0],'-.','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Flange thickness [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Root Flange Breadth [m]','Interpreter', 'latex', 'FontSize',18)
legend('Front Spar','Rear Spar','Interpreter', 'latex', 'FontSize',16)

figure
plot(y,b_front_flange(3,:),'b-','LineWidth',1.5)
hold on
plot(y,b_rear_flange(3,:),'r-','LineWidth',1.5)
plot([0, 1.08, 3.3272],[b_front_flange(3,1),0.01,0.01],'b--','LineWidth',1.5)
plot([0, 1.08, 3.3272],[b_rear_flange(3,1),0.01,0.01],'r--','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Flange Breadth [m]','Interpreter', 'latex', 'FontSize',18)
legend('Front Spar (Ideal)','Rear Spar (Ideal)','Front Spar (Actual)','Rear Spar (Actual)','Interpreter', 'latex', 'FontSize',16)

figure
plot(t_front_flange,m_front_flange,'LineWidth',1.5)
hold on
plot(t_rear_flange,m_rear_flange,'LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Flange thickness [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Flange Mass [kg]','Interpreter', 'latex', 'FontSize',18)
legend('Front Spar','Rear Spar','Interpreter', 'latex', 'FontSize',16)


% figure
% yyaxis left
% plot(t_front_flange,b_front_flange_max,'LineWidth',1.5)
% hold on
% plot(t_rear_flange,b_rear_flange_max,'LineWidth',1.5)
% ylabel('Flange Breadth [m]','Interpreter', 'latex', 'FontSize',18)
% yyaxis right
% plot(t_front_flange,m_front_flange,'LineWidth',1.5)
% hold on
% plot(t_rear_flange,m_rear_flange,'LineWidth',1.5)
% hold off
% grid minor
% box on
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex')
% xlabel('Flange thickness [m]','Interpreter', 'latex', 'FontSize',18)
% ylabel('Flange Mass [kg]','Interpreter', 'latex', 'FontSize',18)
% legend('Front Spar','Rear Spar','Interpreter', 'latex', 'FontSize',16)


%% D-Section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% D-Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_shear_flow = q_front;

set1_1 = [0,8.1; 1.6,9; 2.7,10; 3.35,11; 3.95,12; 4.4,13; 4.8,14; 5.25,15; 6.05,17; 7.1,20; ...
       7.4,21; 8,23; 8.4,24; 10.4,31; 11,33; 11.3,34; 11.8,36; 12.3,38; 12.8,40; 14,45; 14.2,46; ...
       14.9,49; 15.6,52; 15.8,53; 16,54; 16.2,55; 17.1,59];
p1_1 = polyfit(set1_1(:,1)', set1_1(:,2)', 7);

% Model curve ratio = 1.5
set1_15 = [0,6.6; 1.4,7; 2.2,7.5; 2.8,8; 3.55,9; 4.1,10; 4.6,11; 5.1,12; 5.8,13.8; 6.2,15; ...
       6.8,16.5; 7.3,18; 7.9,20; 8.4,21.5; 9,23.5; 9.2,24; 9.5,25; 9.8,26; 10.8,29.5; 11.2,30.9; ...
       11.8,33; 12.6,36; 13.6,40; 14.8,45; 15.7,49; 16.8,54; 17.8,59];
p1_15 = polyfit(set1_15(:,1)', set1_15(:,2)', 7);

% Model curve ratio = 5
set1_5 = [0,4.8; 1,5; 2,5.5; 2.6,6; 3.45,7; 4.4,8.5; 5.4,10.5; 6,12; 6.8,14; 8.95,20; 9.8,23; ...
       10.1,24; 10.4,25; 11.5,29; 12.6,33; 13.6,37; 13.8,38; 14.6,41; 14.8,42; 15.5,45; 16,47; ...
       17.3,53; 18,56; 18.2,56; 18.4,58; 18.8,60];
p1_5 = polyfit(set1_5(:,1)', set1_5(:,2)', 7);

% Plot results
x1 = 0:0.2:20;
y1_1 = polyval(p1_1,x1);
y1_15 = polyval(p1_15,x1);
y1_5 = polyval(p1_5,x1);

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

% Graph 2 (a/b)
% Model curve ratio = 1
set2_1 = [0,8.1; 1.6,9; 2.7,10; 3.35,11; 3.95,12; 4.4,13; 4.8,14; 5.25,15; 6.05,17; 7.1,20; ...
       7.4,21; 8,23; 8.4,24; 10.4,31; 11,33; 11.3,34; 11.8,36; 12.3,38; 12.8,40; 14,45; 14.2,46; ...
       14.9,49; 15.6,52; 15.8,53; 16,54; 16.2,55; 17.1,59];
p2_1 = polyfit(set2_1(:,1)', set2_1(:,2)', 7);

% Model curve ratio = 1.5
set2_15 = [0,6.6; 1.6,6.8; 2.1,7; 2.8,7.5; 3.35,8; 4,9; 4.55,10; 4.8,10.5; 5.4,11.8; 6.3,14; ...
       7.1,16; 7.5,17; 8.8,20.5; 10,24; 10.2,24.5; 10.8,26.5; 11,27; 11.8,29.5; 12.4,31.5; 13,33.5; ...
       14,37; 14.55,39; 15,40.5; 15.7,43; 16.8,47; 17.2,48.5; 17.6,50; 18,51.5; 19.2,56; 19.5,57; ...
       20,58.8];
p2_15 = polyfit(set2_15(:,1)', set2_15(:,2)', 7);

% Model curve ratio = 2
set2_2 = [0,5.8; 2,6.1; 2.4,6.4; 3,7; 3.4,7.5; 3.85,8; 4.4,9; 4.95,10; 5.4,10.9; 6,12.2; ...
       6.6,13.5; 7,14.5; 7.6,15.9; 8.45,18; 9.7,21; 10,22; 11,24.5; 11.2,25; 12.2,28; 13.2,31; ...
       13.8,33; 14.4,35; 15,37; 15.6,39; 16.2,41; 17.2,44.6; 18,47.3; 19,51; 19.6,53; 20,54.5];
p2_2 = polyfit(set2_2(:,1)', set2_2(:,2)', 7);

% Model curve ratio = 3
set2_3 = [0,5.3; 1.4,5.5; 2,5.6; 2.4,5.8; 2.6,6; 3.5,7; 4.6,8.6; 4.8,8.9; 5.4,10; 6,11.2; ...
       6.5,12; 7.4,14; 7.9,15; 8.8,17; 9.2,18; 10,19.8; 10.5,21; 11.3,23; 12.4,26; 13.2,28; ...
       14.3,31; 15,33; 15.4,34; 16.4,37; 16.6,37.5; 17.4,39.8; 17.8,41.1; 18.15,42; 18.6,43.5; ...
       19.1,45; 19.45,46; 19.6,46.5; 19.8,47; 20,47.6];
p2_3 = polyfit(set2_3(:,1)', set2_3(:,2)', 7);

% Model curve ratio = inf
set2_i = [0,4.8; 1.8,5; 2.6,5.5; 3.2,6; 4,7; 4.8,8; 5.4,9; 6,10; 6.6,11; 7.3,12; 7.9,13; ...
       9.15,15; 9.75,16; 10.4,17; 11,18; 12,19.5; 12.6,20.5; 13.2,21.5; 13.6,22.4; 14.2,23; ...
       14.8,24; 16,26; 16.4,26.5; 17,27.5; 18,29; 18.55,30; 19.25,31; 19.6,31.5; 19.85,32; 20,32.3];
p2_i = polyfit(set2_i(:,1)', set2_i(:,2)', 7);

% Plot results
x2 = 0:0.2:20;
y2_1 = polyval(p2_1,x2);
y2_15 = polyval(p2_15,x2);
y2_2 = polyval(p2_2,x2);
y2_3 = polyval(p2_3,x2);
y2_i = polyval(p2_i,x2);

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

% Save polynomial parameters
ratiopoly = [p1_1; p1_15; p1_5; p2_1; p2_15; p2_2; p2_3; p2_i];

tau_tresca = 125*10^6; %%%%%%%%%%%%%%%%%%%%% Update!

figure

for k = 1:5 % iterates the whole dcell calculation 3 times to to get converged results
    n_pseudo_ribs = 0:15;
    a_pseudo_ribs = s./((n_pseudo_ribs + 1)*cosd(sweep_LE));
    t_range = [0.001,0.0015,0.002,0.0025,0.003];
    t_dcell = t_range(k)*ones(1,length(n_pseudo_ribs));

    data = [0.14645,  0.03994;0.12408,  0.03795;0.10332,  0.03564;0.08427,  0.03305;0.06699,  0.03023;0.05156,  0.02720;0.03806,  0.02395;0.02653,  0.02039;0.01704,  0.01646;0.00961,  0.01214;0.00428,  0.00767;0.00107,  0.00349;0.0,      0.0];

    l = 0;
    for j = 1:length(y)
        for i = 1:12
            l(i,j) = sqrt((c(j)*data(i+1,1) - c(j)*data(i,1))^2 + (c(j)*data(i+1,2) - c(j)*data(i,2))^2);
            A(i,j) = trapz([c(j)*data(i,1),c(j)*data(i+1,1)],[c(j)*data(i,2),c(j)*data(i+1,2)]);
        end
        l_curved(j) = (2*sum(l(:,j)) + 2*(0.15-0.14645)*c(j))/2;  
        A_pseudo_rib(j) = -2*sum(A(:,j));
    end

    Radius = (((x_FS.*c).^2)./(y_avg/2)+((y_avg/2).^2)./(x_FS.*c))/2;   % 0.9684

    for i = 1:length(n_pseudo_ribs)
        a_b(i,:) = a_pseudo_ribs(i)./l_curved;
        b_a(i,:) = l_curved/a_pseudo_ribs(i);
    end

    for i = 1:length(n_pseudo_ribs)
        for j = 1:length(y)
            if a_b(i,j) >= 1
                l_sqrt_Rt(i,j) = l_curved(j)/sqrt(Radius(j)*t_dcell(i)); %b_sqrt_Rt
            else
                l_sqrt_Rt(i,j) = a_pseudo_ribs(i)/sqrt(Radius(j)*t_dcell(i)); %a_sqrt_Rt
            end
        end
    end
    for i = 1:length(n_pseudo_ribs)
        Ks_dcell(i,:) = DSectionRatioFunc(b_a(i,:),l_sqrt_Rt(i,:),ratiopoly);
        tau_crit(i,:) = Ks_dcell(i,:).*E.*t_dcell(i).^3./(l_curved.^2);
        tau_crit_min(i) = min(tau_crit(i,:));
%         t_panel(i,:) = sqrt(tau_tresca.*l_curved.^2./(Ks_dcell(i,:).*E));
%         t_opt(i) = max(t_panel(i,:));
    end

    plot(n_pseudo_ribs,tau_crit_min,'LineWidth',1.5)
    hold on
end
plot([0 15],[max_shear_flow(1),max_shear_flow(1)],'k--','LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Number of Pseudo Ribs','Interpreter', 'latex', 'FontSize',18)
ylabel('Minimum Shear Flow [N/m]','Interpreter', 'latex', 'FontSize',18)
legend('$t = 1\,mm$','$t = 1.5\,mm$','$t = 2\,mm$','$t = 2.5\,mm$','$t = 3\,mm$','Max Applied Shear Flow','Interpreter', 'latex', 'FontSize',18)

m_dcell = rho*t_range*a_pseudo_ribs(1)*(l_curved(1)+l_curved(end))/2;
t_pseudo_rib = 0.001;

for i = 1:length(n_pseudo_ribs)
    m_pseudo_ribs(i) = 0;
    if n_pseudo_ribs(i) > 0
        for j = 1:n_pseudo_ribs(i)
            m_pseudo_ribs(i) = m_pseudo_ribs(i) + rho*A_pseudo_rib(floor((length(y)*a_pseudo_ribs(i)*j)/(s/cosd(sweep_LE))))*t_pseudo_rib*n_pseudo_ribs(i);
        end
    end
end

for i = 1:length(t_range)
    for j = 1:length(n_pseudo_ribs)
        m_total_dcell(i,j) = m_dcell(i) + m_pseudo_ribs(j);
    end
end

figure
plot(n_pseudo_ribs,m_total_dcell(1,:),'LineWidth',1.5)
hold on
plot(n_pseudo_ribs,m_total_dcell(2,:),'LineWidth',1.5)
plot(n_pseudo_ribs,m_total_dcell(3,:),'LineWidth',1.5)
plot(n_pseudo_ribs,m_total_dcell(4,:),'LineWidth',1.5)
plot(n_pseudo_ribs,m_total_dcell(5,:),'LineWidth',1.5)
hold off
grid on
box on
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Number of Pseudo Ribs','Interpreter', 'latex', 'FontSize',18)
ylabel('Overall D-Cell Mass [kg]','Interpreter', 'latex', 'FontSize',18)
legend('$t = 1\,mm$','$t = 2\,mm$','$t = 3\,mm$','Interpreter', 'latex', 'FontSize',18)


%% Wing Configuration

%%%%%%%%%%%%%%%%%%%%%%%%% Plot Wing Configuration %%%%%%%%%%%%%%%%%%%%%%%%

n = 25; % number of stringers
n_pseudo = 0;

figure
p1 = plot([0,s,s,0,0],[c_root,c_root-s*tand(sweep_LE),c_root-s*tand(sweep_LE)-c_tip,0,c_root],'k-','LineWidth',1.5);
hold on 
plot([0,s],[c_root*(1-x_RS),c_tip*(1-x_RS)+c_root-s*tand(sweep_LE)-c_tip],'b-','LineWidth',1.5)
z_LE = @(x) ((c_root-(c_root-s*tand(sweep_LE)))/(0-s))*x + c_root;
z = @(x) ((c_root*(1-x_FS)-(c_tip*(1-x_FS)+c_root-s*tand(sweep_LE)-c_tip))/(0-s))*x + c_root*(1-x_FS);
for i = 1:n
    hoe = c_tip*(1-x_RS)+c_root-s*tand(sweep_LE)-c_tip+i*b_stringers(n+1);
    if hoe <= c_tip*(1-x_FS)+c_root-s*tand(sweep_LE)-c_tip
        x_st(i,:) = [0,s];
        y_st(i,:) = [c_root*(1-x_RS)+i*b_stringers(n+1),hoe];
        p3 = plot(x_st(i,:),y_st(i,:),'color',[0, 0.40, 0],'LineWidth',1);
    else
        z_stringer = @(x) ((c_root*(1-x_RS)+i*b_stringers(n+1)-hoe)/(0-s))*x + c_root*(1-x_RS)+i*b_stringers(n+1);
        Int = fzero(@(x) z(x)-z_stringer(x), 1);
        plot([0,Int],[c_root*(1-x_RS)+i*b_stringers(n+1),z(Int)],'color',[0, 0.40, 0],'LineWidth',1)
    end
end
for j = 1:size(RibLocations,2)-2
    if RibLocations(i,j) ~= 0
        p4 = plot([0+y(RibLocations(n+1,j)),0+y(RibLocations(n+1,j))],[z(y(RibLocations(n+1,j))),z(y(RibLocations(n+1,j)))-l_box(RibLocations(n+1,j))],'r-','LineWidth',1.5);
        
    end
end
if n_pseudo > 0
    for k = 1:n_pseudo
        spacing = floor(length(y)/(n_pseudo+1));
        p5 = plot([0+y(spacing*k),0+y(spacing*k)],[z_LE(y(spacing*k)),z(y(spacing*k))],'Color',[0.75, 0.75, 0],'LineWidth',1.5);
    end
end
p2 = plot([0,s],[c_root*(1-x_FS),c_tip*(1-x_FS)+c_root-s*tand(sweep_LE)-c_tip],'b-','LineWidth',1.5);
p6 = fill([0,0,s*0.9,s*0.9,0],[0,c_root*0.32,c_root*0.32-s*0.9*tand(sweep_hinge),c_root*0.32-s*0.9*tand(sweep_hinge)-1.01,0],[0.827,0.827,0.827]); %0 0.9 0.5
p1 = plot([0,0,s*0.9,s*0.9,0],[0,c_root*0.32,c_root*0.32-s*0.9*tand(sweep_hinge),c_root*0.32-s*0.9*tand(sweep_hinge)-1.01,0],'k-','LineWidth',1);
hold off
box on 
ax=gca;ax.LineWidth=1;
set(gca,'FontSize',15)
set(gca,'TickLabelInterpreter','latex')
xlabel('Spanwise Position [m]','Interpreter', 'latex', 'FontSize',18)
ylabel('Chordwise Position [m]','Interpreter', 'latex', 'FontSize',18)
legend([p2 p3 p4 p6],{'Spars','Stringers','Ribs','Rudder'},'Interpreter', 'latex', 'FontSize',16)


% figure
% hold on
% pos = 0;
% for i = 1:5
%     surf([0 s;0 s],[0+pos 0+pos;0+pos 0+pos],[0 0;h h],'LineWidth',2)
%     surf([0 s;0 s],[0+pos 0+pos;0+pos+d_h_CT*h 0+pos+d_h_CT*h],[0 0;0 0],'LineWidth',2)
%     surf([0 s;0 s],[0+pos 0+pos;0+pos-d_h_CT*h 0+pos-d_h_CT*h],[h h;h h],'LineWidth',2)
%     pos = pos + b_stringers(27);
% end
% surf([0 s;0 s],[0-b_stringers(27) 0-b_stringers(27);0+pos 0+pos],[0 0;0 0],'LineWidth',1)
% hold off
% grid on
% box on
% ax=gca;ax.LineWidth=1;
% set(gca,'FontSize',15)
% set(gca,'TickLabelInterpreter','latex')
% ylim([-1 1])
% zlim([-1 1])
% % xlabel('Number of Stringers','Interpreter', 'latex', 'FontSize',18)
% % ylabel('Stringer Height [\%]','Interpreter', 'latex', 'FontSize',18)
% % zlabel('Mass [kg]','Interpreter', 'latex', 'FontSize',18)
% view(40,40);
