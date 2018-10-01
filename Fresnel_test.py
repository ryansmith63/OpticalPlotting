from TSTR_fit_new import F

T_330=1-F(0.,1.0,1.48,0.5)*4
T_190=1-F(0.,1.0,1.55,0.5)*4
T_170=1-F(0.,1.0,1.62,0.5)*4
print("Transmission, 330 nm: {0:g}, 190 nm: {1:g}; ratio: {2:g}".format(T_330,T_190,T_330/T_190))