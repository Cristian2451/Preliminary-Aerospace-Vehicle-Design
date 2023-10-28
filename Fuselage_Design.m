%Fuselage Design
clear
clc

%Plane Values
d = 2.786;

%Material Values
sigma_y = 290e6; %Check this value using Al 2024-T3


%Pressure loads
P_int = 0.7e5;
P_ext = 18753.93;
SafetyFactor = 1.5;

P_load = SafetyFactor*(P_int - P_ext);

%Defined by longitudinal pressure stress
t_skin_min = (d * P_load)/(4 * sigma_y);

%Matching strains
t_bulk_min = 2.5 * t_skin_min;


%Stringer Info
n_s = 18;              %First guess at the number of stringers
alpha = 2 * pi / n_s;

      %Co-ords of the stringers
z = 0.5 * d * cos(alpha * [1:n_s]);
y = 0.5 * d * sin(alpha * [1:n_s]);

plot(y,z)
axis equal
b = 0.5 * d * alpha; %arc length between each stringer

%Idealising as a Boom Shear Panel to get a first pass
BM_max = -1019935; %Max bending moment
t = t_skin_min;

A_s_initial = abs(BM_max) / (0.5 * d * n_s * sigma_y); %scalar value, constant for all

A_s = A_s_initial * [1:n_s];

        %These are used for matrix boom stuff
z1 = 0.5 * d * cos(alpha * [[2:n_s],1]); %Shift of z, 1 to the left (zi+1)
z2 = 0.5 * d * cos(alpha * [n_s,[1:n_s-1]]); %Shift of z 1 to the right (zi-1)

B = (b * t / 6)*(4 + (z1+z2)./z) + A_s;     %Boom area given by boom modelling equation

Iyy = sum(z.^2 .* B);