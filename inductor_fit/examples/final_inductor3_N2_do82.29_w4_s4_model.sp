* WIDEBAND RFIC CENTER-TAPPED INDUCTOR MODEL WITH SUBSTRATE NETWORK
* Created using inductor_fit.py from final_inductor3_N2_do82.29_w4_s4.s3p
* Fit band 0.010 to 60.000 GHz, 2 coil segment(s) per half coil
* k is the effective coupling including the center tap lead inductance
* Total element values:
*   Rs_h1 = 0.721296
*   Ls_h1 = 1.47676e-10
*   Rskin1_h1 = 0.816092
*   Lskin1_h1 = 1.12458e-11
*   Rskin2_h1 = 13.7486
*   Lskin2_h1 = 1.35463e-11
*   Rs_h2 = 0.998252
*   Ls_h2 = 1.47676e-10
*   Rskin1_h2 = 1.18808
*   Lskin1_h2 = 1.12458e-11
*   Rskin2_h2 = 19.8508
*   Lskin2_h2 = 1.35463e-11
*   k = 0.423279
*   Cs = 4.8286e-15
*   Rct = 0.0957937
*   Cox1 = 6.47091e-15
*   Rsi1 = 524.74
*   Csi1 = 5.59017e-14
*   Cox2 = 7.17834e-15
*   Rsi2 = 561.036
*   Csi2 = 5.38782e-14
*   Coxct = 1.90788e-14
*   Rsict = 8364.33
*   Csict = 2.91126e-15
*
.SUBCKT ct_inductor_model p1 p2 ct
* half coil 1, segment 1
Rs_h1_1 p1 a_h1_1 0.360648
Ls_h1_1 a_h1_1 b_h1_1 7.3838e-11
Rskin1_h1_1 b_h1_1 d_h1_1 0.408046
Lskin1_h1_1 b_h1_1 d_h1_1 5.62289e-12
Rskin2_h1_1 d_h1_1 c1 6.87429
Lskin2_h1_1 d_h1_1 c1 6.77314e-12
* half coil 1, segment 2
Rs_h1_2 c1 a_h1_2 0.360648
Ls_h1_2 a_h1_2 b_h1_2 7.3838e-11
Rskin1_h1_2 b_h1_2 d_h1_2 0.408046
Lskin1_h1_2 b_h1_2 d_h1_2 5.62289e-12
Rskin2_h1_2 d_h1_2 m 6.87429
Lskin2_h1_2 d_h1_2 m 6.77314e-12
* half coil 2, segment 1
Rs_h2_1 m a_h2_1 0.499126
Ls_h2_1 a_h2_1 b_h2_1 7.3838e-11
Rskin1_h2_1 b_h2_1 d_h2_1 0.594039
Lskin1_h2_1 b_h2_1 d_h2_1 5.62289e-12
Rskin2_h2_1 d_h2_1 c3 9.9254
Lskin2_h2_1 d_h2_1 c3 6.77314e-12
* half coil 2, segment 2
Rs_h2_2 c3 a_h2_2 0.499126
Ls_h2_2 a_h2_2 b_h2_2 7.3838e-11
Rskin1_h2_2 b_h2_2 d_h2_2 0.594039
Lskin1_h2_2 b_h2_2 d_h2_2 5.62289e-12
Rskin2_h2_2 d_h2_2 p2 9.9254
Lskin2_h2_2 d_h2_2 p2 6.77314e-12
* magnetic coupling of the half coils (dot at the first node, aiding for differential current)
K_h1_1_h2_1 Ls_h1_1 Ls_h2_1 0.21164
K_h1_1_h2_2 Ls_h1_1 Ls_h2_2 0.21164
K_h1_2_h2_1 Ls_h1_2 Ls_h2_1 0.21164
K_h1_2_h2_2 Ls_h1_2 Ls_h2_2 0.21164
Cs p1 p2 4.8286e-15
Rct m ct 0.0957937
* substrate network at p1
Cox_n0 p1 s0 3.23546e-15
Rsi_n0 s0 0 1049.48
Csi_n0 s0 0 2.79509e-14
* substrate network at c1
Cox_n1 c1 s1 8.00517e-15
Rsi_n1 s1 0 1017.56
Csi_n1 s1 0 2.86787e-14
* substrate network at m
Cox_n2 m s2 9.53942e-15
Rsi_n2 s2 0 16728.7
Csi_n2 s2 0 1.45563e-15
* substrate network at c3
Cox_n3 c3 s3 8.35888e-15
Rsi_n3 s3 0 1085.66
Csi_n3 s3 0 2.76669e-14
* substrate network at p2
Cox_n4 p2 s4 3.58917e-15
Rsi_n4 s4 0 1122.07
Csi_n4 s4 0 2.69391e-14
.ENDS
