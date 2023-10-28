function [mass] = lowerSurface(N,h,ts_t, yield,BM,L, ribLoc,M,Loc)
format long;
E = 73.8e9; %72.4e9;
trig = false;

root_chord = 4.027;
box_width_root = (0.6-0.15)*root_chord;
web_flange = 0.3; 
pannel_width = box_width_root / N;
h_b = h / pannel_width;
h_root = h;
d_root = (web_flange * h_root);

% Check that K > 3.62 at root
Kroot = K_interp(h_b, ts_t,M);

% Wing box modification to reduce NC
% Fuselage diameter = 2.786
modL = 2.786/2+ 0.2;
boxW_root= 1.9;
boxH_root = 0.4;

boxW_mod = (4.027 - 0.1713.*modL) * (0.6-0.15);
boxH_mod = (4.027 - 0.1713.*modL) * (0.0989+0.0879)/2;
mMod = (boxH_root - boxH_mod) / -modL;
mModW = (boxW_root - boxW_mod) / -modL;

y = 0;

% ts is determined by max stress location, which now at interface of
% changing box structures
box_heightJoin = boxH_root + mMod * y;
box_widthJoin = boxW_root + mModW * y;
NC = interp1(Loc, BM, modL) / (box_widthJoin * box_heightJoin);
t2max = (NC * (pannel_width^2) / (3.62 * E)) ^ (1/3);
ts_set = ts_t * t2max;
i=1;

for i = 1: length(ribLoc)
    chord = (4.027 - 0.1713*y);
    if y < modL
        box_height = boxH_root + mMod * y;
        box_width = boxW_root + mModW * y;
    else
        box_height = chord*(0.0989+0.0879)/2;
        box_width = (0.6-0.15)*chord;
    end

    NC = interp1(Loc, BM, y) / (box_width * box_height);
    % t2 is set using local buckling of a pannel
    t2 = abs((NC * (pannel_width^2) / (3.62 * E)) ^ (1/3)); 

    if t2 < 0.001
        t2_real = 0.001;
    else
        t2_real = t2;
    end
    
    if i == 1
        % If at the root, calculate stringer thickness
        if t2max > t2
            ts_root = ts_set;
        else
            ts_root = ts_t * t2;
        end
        As_root = (h_root + 2*d_root)*ts_root;
    end

    T_eq = t2_real + As_root/pannel_width;
    T_eqList(i) = T_eq;

    actual_stress = NC / T_eq;
    if actual_stress > yield
        mass = 1;
        trig = true;
        break
    end

    K = K_interp(h_b, ts_t,M);
    sigma_pannels = K * E * (t2_real/ pannel_width)^2;

    I_A = h^2 * h_b * ts_t * (0.633+ 0.37* h_b* ts_t) / (1.6 *h_b * ts_t + 1)^2;
    F = ((K * (t2/pannel_width)^2 * pi^2 * I_A)^(1/4) ) / (t2* (1 + 1.6*h_b * ts_t))^(1/2);

    sigma_stringers = F * sqrt(abs(NC * E / L(i)));

    if sigma_pannels < actual_stress
        mass = 2; % Pannels buckle below applied stress
        trig = true;
        break
    elseif sigma_stringers < actual_stress
        mass = 3; % Stringers buckle below applied stress
        trig = true;
        break
    end

    %i = i +1;
    y = ribLoc(i);
end

if Kroot < 3.62
    mass = 4; 
elseif trig == false
    ribLoc2 = [ribLoc, 14.1];
%     disp(length(ribLoc2))
%     disp(length(N))
%     disp(length(T_eqList))
    pannelVol = PannelsVol(N,ribLoc2,T_eqList);
    mass = pannelVol * 2770 ;%2800;
    %As_bt_final = As_bt;
end

end