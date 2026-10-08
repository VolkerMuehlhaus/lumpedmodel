* WIDEBAND RFIC CENTER-TAPPED INDUCTOR MODEL WITH SUBSTRATE NETWORK
* Created using inductor_fit.py from final_inductor3_N2_do82.29_w4_s4.s3p
* Fit band 0.010 to 60.000 GHz, 2 coil segment(s) per half coil
* k is the effective coupling including the center tap lead inductance
* Total element values:
*   Rs_h1 = 0.719718
*   Ls_h1 = 1.51808e-10
*   Rskin1_h1 = 0.713249
*   Lskin1_h1 = 1.09128e-11
*   Rskin2_h1 = 6.10196
*   Lskin2_h1 = 1.00813e-11
*   Rs_h2 = 0.998336
*   Ls_h2 = 1.51808e-10
*   Rskin1_h2 = 1.08414
*   Lskin1_h2 = 1.09128e-11
*   Rskin2_h2 = 9.42995
*   Lskin2_h2 = 1.00813e-11
*   k = 0.411631
*   Cs = 4.84257e-15
*   Rct = 0.0956473
*   Cox1 = 6.48035e-15
*   Rsi1 = 526.721
*   Csi1 = 5.57751e-14
*   Cox2 = 7.18597e-15
*   Rsi2 = 562.284
*   Csi2 = 5.37876e-14
*   Coxct = 1.90461e-14
*   Rsict = 8348.82
*   Csict = 2.92584e-15
*
.SUBCKT ct_inductor_model p1 p2 ct
* half coil 1, segment 1
Rs_h1_1 p1 a_h1_1 0.359859
Ls_h1_1 a_h1_1 b_h1_1 7.59041e-11
Rskin1_h1_1 b_h1_1 d_h1_1 0.356625
Lskin1_h1_1 b_h1_1 d_h1_1 5.45641e-12
Rskin2_h1_1 d_h1_1 c1 3.05098
Lskin2_h1_1 d_h1_1 c1 5.04066e-12
* half coil 1, segment 2
Rs_h1_2 c1 a_h1_2 0.359859
Ls_h1_2 a_h1_2 b_h1_2 7.59041e-11
Rskin1_h1_2 b_h1_2 d_h1_2 0.356625
Lskin1_h1_2 b_h1_2 d_h1_2 5.45641e-12
Rskin2_h1_2 d_h1_2 m 3.05098
Lskin2_h1_2 d_h1_2 m 5.04066e-12
* half coil 2, segment 1
Rs_h2_1 m a_h2_1 0.499168
Ls_h2_1 a_h2_1 b_h2_1 7.59041e-11
Rskin1_h2_1 b_h2_1 d_h2_1 0.542072
Lskin1_h2_1 b_h2_1 d_h2_1 5.45641e-12
Rskin2_h2_1 d_h2_1 c3 4.71498
Lskin2_h2_1 d_h2_1 c3 5.04066e-12
* half coil 2, segment 2
Rs_h2_2 c3 a_h2_2 0.499168
Ls_h2_2 a_h2_2 b_h2_2 7.59041e-11
Rskin1_h2_2 b_h2_2 d_h2_2 0.542072
Lskin1_h2_2 b_h2_2 d_h2_2 5.45641e-12
Rskin2_h2_2 d_h2_2 p2 4.71498
Lskin2_h2_2 d_h2_2 p2 5.04066e-12
* magnetic coupling of the half coils (dot at the first node, aiding for differential current)
K_h1_1_h2_1 Ls_h1_1 Ls_h2_1 0.205816
K_h1_1_h2_2 Ls_h1_1 Ls_h2_2 0.205816
K_h1_2_h2_1 Ls_h1_2 Ls_h2_1 0.205816
K_h1_2_h2_2 Ls_h1_2 Ls_h2_2 0.205816
Cs p1 p2 4.84257e-15
Rct m ct 0.0956473
* substrate network at p1
Cox_n0 p1 s0 3.24017e-15
Rsi_n0 s0 0 1053.44
Csi_n0 s0 0 2.78876e-14
* substrate network at c1
Cox_n1 c1 s1 8.0017e-15
Rsi_n1 s1 0 1021.23
Csi_n1 s1 0 2.8619e-14
* substrate network at m
Cox_n2 m s2 9.52305e-15
Rsi_n2 s2 0 16697.6
Csi_n2 s2 0 1.46292e-15
* substrate network at c3
Cox_n3 c3 s3 8.35451e-15
Rsi_n3 s3 0 1087.93
Csi_n3 s3 0 2.76252e-14
* substrate network at p2
Cox_n4 p2 s4 3.59299e-15
Rsi_n4 s4 0 1124.57
Csi_n4 s4 0 2.68938e-14
.ENDS
