
moment = 742.39*span_loc^3 -33086*span_loc^2 +494596.53*span_loc -2473450.12;

if span_loc<=4
    shear=-975.16*span_loc^2 +45163*span_loc-459499;
else
    shear=-53.11*span_loc^3-26.682*span_loc^2+38704*span_loc-389754;
end

torque=[-4970.901227
-4857.927095
-4739.521194
-4616.032632
-4487.791686
-4355.111281
-4218.288146
-4077.603665
141971.8028
-3785.702343
-3634.974634
-3481.363018
-3325.072594
-3166.28988
-3005.179867
-2841.881696
-2676.502335
-2509.10731
-2339.706923
-2168.235321
-1994.517758
-1818.217375
-1638.744333
-1455.090276
-1265.499553
-1066.731454
-852.0628655
-594.131857
-0.560303572];